function dy = finite_diff_central(y, t_or_dt)
%FINITE_DIFF_CENTRAL  Numerical derivative using the central difference method.
%
%   dy = finite_diff_central(y)
%       Assumes uniformly sampled data with unit spacing (dt = 1).
%
%   dy = finite_diff_central(y, dt)
%       Assumes uniform sampling with constant time step dt (scalar).
%       Example: dy = diff_central(f, 0.1)  % sampling interval = 0.1 s
%
%   dy = finite_diff_central(y, t)
%       Uses the explicit time vector t (size N×1) to compute derivatives
%       with potentially nonuniform sampling. The vector t must be strictly
%       increasing and have the same number of elements as the rows of y.
%
% INPUTS
%   y        [N×D] matrix of function samples (columns are independent series)
%   t_or_dt  either:
%             - scalar  → uniform step dt
%             - vector  → time (or independent variable) samples [N×1]
%
% OUTPUT
%   dy       [N×D] matrix of first derivative estimates.
%
% METHOD
%   Central differences for interior points (second-order accuracy, O(dt²)).
%   Forward and backward differences (first-order) are used at the endpoints.
%
% NOTE
%   This function can be used to estimate the derivative of any function f(x)
%   from its sampled values at a finite, discrete set of points in the
%   independent variable (uniform or nonuniform sampling).
%
% v2.0 - Francesco Fiorina - optimized for loops

if nargin < 2
    t_or_dt = 1;
end

[N,D] = size(y);
dy = nan(N,D);

% ---------- UNIFORM SAMPLING ----------
if isscalar(t_or_dt)

    dt = t_or_dt;

    if dt <= 0
        error('dt must be positive.')
    end

    % central difference
    dy(2:N-1,:) = (y(3:N,:) - y(1:N-2,:)) / (2*dt);

    % edges
    dy(1,:) = (y(2,:) - y(1,:)) / dt;
    dy(N,:) = (y(N,:) - y(N-1,:)) / dt;

% ---------- NON-UNIFORM SAMPLING ----------
else

    t = t_or_dt(:);

    if numel(t) ~= N
        error('Time vector length mismatch')
    end

    if any(diff(t)<=0)
        error('Time must be strictly increasing')
    end

    % central difference
    dt = t(3:N) - t(1:N-2);
    dy(2:N-1,:) = (y(3:N,:) - y(1:N-2,:)) ./ dt;

    % edges
    dy(1,:) = (y(2,:) - y(1,:)) / (t(2)-t(1));
    dy(N,:) = (y(N,:) - y(N-1,:)) / (t(N)-t(N-1));

end



%% - Old version - deprecated
% if nargin < 2
%     t_or_dt = 1;   % default step = 1 
% end
% 
% [N, D] = size(y);
% dy = nan(N, D);
% 
% % --- Scalar dt ---
% if isscalar(t_or_dt)
%     
%     dt = t_or_dt;
%     if dt <= 0, error('dt must be positive.'); end
%     % Loop over columns (independent series) --> ROOM FOR OPTIMIZATION!
%     for j = 1:D
%         f = y(:,j);
%         % central differences
%         dy(2:N-1,j) = (f(3:N) - f(1:N-2)) / (2*dt);
%         % forward/backward at edges
%         dy(1,j) = (f(2) - f(1)) / dt;
%         dy(N,j) = (f(N) - f(N-1)) / dt;
%     end
% 
% % --- Vector t ---    
% else
%     
%     t = t_or_dt(:);
%     if numel(t) ~= N
%         error('The time vector t must have %d elements.', N);
%     end
%     if any(~isfinite(t)) || any(diff(t) <= 0)
%         error('The time vector must be finite and strictly increasing.');
%     end
%     % Loop over columns (independent series) --> ROOM FOR OPTIMIZATION!
%     for j = 1:D
%         f = y(:,j);
%         % central differences
%         for i = 2:N-1
%             dy(i,j) = (f(i+1)-f(i-1)) / (t(i+1)-t(i-1));
%         end
%         % forward/backward at edges
%         dy(1,j) = (f(2)-f(1)) / (t(2)-t(1));
%         dy(N,j) = (f(N)-f(N-1)) / (t(N)-t(N-1));
%     end
% 
% end
% 
% end