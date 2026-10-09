function Y = movstd_nancheck(X, K, nan_threshold, varargin)
%MOVSTD_NANCHECK  Moving standard deviation with NaN validity threshold.
%
%   Y = MOVSTD_NANCHECK(X, K, NAN_THRESHOLD) computes a centered moving
%   standard deviation along vector (or N-D array) X using a window of
%   length K, with an additional per-window NaN check.
%
%   NAN_THRESHOLD is a value in (0, 1] that represents the *maximum allowed
%   fraction* of NaN values inside a window:
%
%       - If  (number of NaNs in window) / (window size)  >  NAN_THRESHOLD
%         the output sample Y(i) is set to NaN, regardless of the data.
%
%       - Otherwise the standard deviation is computed using the non-NaN
%         elements of the window (equivalent to 'omitnan').
%
%   This lets you write, for example:
%       NAN_THRESHOLD = 0.5   ->  a window must have at least 50 % valid
%                                 samples; otherwise the result is NaN.
%       NAN_THRESHOLD = 0.0   ->  any single NaN contaminates the window
%                                 (same behaviour as 'includenan').
%       NAN_THRESHOLD = 1.0   ->  NaN values are always omitted; the result
%                                 is NaN only when *all* window values are
%                                 NaN (same behaviour as 'omitnan').
%
%   Y = MOVSTD_NANCHECK(X, K, NAN_THRESHOLD) uses the default centred
%   window. K must be a positive integer scalar or a two-element vector
%   [NB NF] of non-negative integers (same semantics as MOVSTD).
%
%   Y = MOVSTD_NANCHECK(X, K, NAN_THRESHOLD, NRM) additionally specifies
%   the normalisation:
%       0  - (default) normalise by (n-1), where n is the number of valid
%            samples used in each window.
%       1  - normalise by n.
%
%   Y = MOVSTD_NANCHECK(X, K, NAN_THRESHOLD, NRM, DIM) operates along
%   dimension DIM of X.
%
%   Notes
%   -----
%   * If the centre sample X(i) itself is NaN, Y(i) is set to NaN
%     immediately and the window computation is skipped entirely.
%   * Endpoint handling always uses the 'shrink' strategy (the window is
%     reduced to fit at the borders), consistent with the MOVSTD default.
%   * NAN_THRESHOLD is applied to the *effective* window size at each
%     position (i.e. after shrinking at the endpoints).
%   * The function works for vectors and N-D arrays.
%
%   Examples
%   --------
%   % Basic usage – allow up to 40 % NaNs per window
%   x = [4 8 NaN -1 -2 NaN -1 3 4 5];
%   y = movstd_nancheck(x, 5, 0.4);
%
%   % Stricter threshold – any NaN invalidates the window
%   y_strict = movstd_nancheck(x, 5, 0.0);
%
%   % Permissive – always omit NaNs (equivalent to movstd 'omitnan')
%   y_permissive = movstd_nancheck(x, 5, 1.0);
%
%   % Two-sided asymmetric window [NB NF] = [3 1], normalise by n
%   y_asym = movstd_nancheck(x, [3 1], 0.5, 1);
%
%   See also MOVSTD, STD, MOVMEAN, MOVVAR
% -------------------------------------------------------------------------
%   AUTHOR:
%       Francesco Fiorina, NavSAS - DET - Politecnico di Torino, 04/03/2026
%
%   VERSION:
%       1.0
%
%   COPYRIGHT:
%       NavSAS - DET - Politecnico di Torino, 2026. All rights reserved.
% -------------------------------------------------------------------------

    % ------------------------------------------------------------------
    % 1.  Input parsing
    % ------------------------------------------------------------------
    narginchk(3, 5);

    % NRM (normalisation flag) – default 0
    nrm = 0;
    if nargin >= 4 && ~isempty(varargin{1})
        nrm = varargin{1};
        validateattributes(nrm, {'numeric','logical'}, {'scalar','binary'}, ...
            'movstd_nancheck', 'NRM');
    end

    % DIM – default: first non-singleton dimension
    sz = size(X);
    dim = find(sz ~= 1, 1, 'first');
    if isempty(dim), dim = 1; end          % X is a scalar
    if nargin >= 5 && ~isempty(varargin{2})
        dim = varargin{2};
        validateattributes(dim, {'numeric'}, ...
            {'scalar','integer','positive','<=',ndims(X)}, ...
            'movstd_nancheck', 'DIM');
    end

    % NAN_THRESHOLD validation
    validateattributes(nan_threshold, {'numeric'}, ...
        {'scalar','real','>=',0,'<=',1}, 'movstd_nancheck', 'NAN_THRESHOLD');

    % Window specification  -----------------------------------------------
    if isscalar(K)
        validateattributes(K, {'numeric'}, {'scalar','integer','positive'}, ...
            'movstd_nancheck', 'K');
        if mod(K,2) == 1          % odd  -> symmetric
            nb = (K-1)/2;
            nf = nb;
        else                      % even -> one extra past element (MATLAB convention)
            nb = K/2;
            nf = K/2 - 1;
        end
    elseif numel(K) == 2
        nb = K(1);  nf = K(2);
        validateattributes(nb, {'numeric'}, {'scalar','integer','nonnegative'}, ...
            'movstd_nancheck', 'K(1)=NB');
        validateattributes(nf, {'numeric'}, {'scalar','integer','nonnegative'}, ...
            'movstd_nancheck', 'K(2)=NF');
        K  = nb + nf + 1;         % total nominal window length
    else
        error('movstd_nancheck:badK', ...
            'K must be a positive integer scalar or a two-element [NB NF] vector.');
    end

    % ------------------------------------------------------------------
    % 2.  Work along the chosen dimension by permuting so dim -> dim 1
    % ------------------------------------------------------------------
    perm  = [dim, setdiff(1:ndims(X), dim)];   % permutation vector
    Xp    = permute(X, perm);                  % dim of interest is now dim 1
    szp   = size(Xp);
    N     = szp(1);                            % length along working dimension
    nOther = prod(szp(2:end));                 % all other elements

    Xmat  = reshape(Xp, N, nOther);            % 2-D: [N x rest]
    Ymat  = NaN(N, nOther);                    % pre-allocate output

    % ------------------------------------------------------------------
    % 3.  Slide the window
    % ------------------------------------------------------------------
    for i = 1:N
        % --- Fast-path: skip entirely if the centre sample is NaN --------
        % Ymat is pre-allocated to NaN, so skipped positions are already set.
        centre_nan = isnan(Xmat(i, :));   % [1 x nOther] logical
        if all(centre_nan)
            continue                       % every column is NaN at centre -> skip
        end

        % Effective window indices (shrink at endpoints)
        i_start = max(1,   i - nb);
        i_end   = min(N,   i + nf);
        win_len = i_end - i_start + 1;   % effective (possibly shrunk) window

        % Extract window for every "column" (remaining dimensions)
        win = Xmat(i_start:i_end, :);    % [win_len x nOther]

        % Count NaNs per column
        n_nan    = sum(isnan(win), 1);   % [1 x nOther]
        frac_nan = n_nan ./ win_len;     % fraction of NaNs

        % valid_cols is a logical [1 x nOther] mask (scalar for 1-D input).
        % No find() on a scalar flag – we iterate only over truly valid columns.
        valid_cols = ~(frac_nan > nan_threshold | centre_nan);
 
         for c = find(valid_cols)
             vals = win(~isnan(win(:, c)), c);   % non-NaN elements of this column
             n    = numel(vals);
             if n == 0
                 % All elements were NaN – safety net, should not reach here
                 % after threshold check, but kept for robustness.
                 Ymat(i, c) = NaN;
             elseif n == 1
                 % Single valid sample: std is 0 regardless of normalisation.
                 % This matches movstd(...,'omitnan') behaviour at shrunk
                 % endpoints (e.g. first/last element of the vector).
                 Ymat(i, c) = 0;
             else
                 Ymat(i, c) = std(vals, nrm);    % nrm=0 -> /n-1 ; nrm=1 -> /n
             end
         end
         % Columns where valid_cols is false remain NaN (pre-allocated)
    end

    % ------------------------------------------------------------------
    % 4.  Reshape and undo permutation
    % ------------------------------------------------------------------
    Yp = reshape(Ymat, szp);
    % Inverse permutation
    iperm(perm) = 1:numel(perm);
    Y = ipermute(Yp, perm);

end