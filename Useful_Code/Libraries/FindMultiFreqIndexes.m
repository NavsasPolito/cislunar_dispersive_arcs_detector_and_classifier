function [prnIdx,obsSolutions] = FindMultiFreqIndexes(obsSolutions)
% This function findes the index of each PRN on each band
%
% Written by Simone Zocca

%% --- Keep a list of PRNs present for each of the three bands
prnLen = 400;
prnFullIdx = NaN(3,prnLen);
prnIdx = NaN(2,prnLen);

% Offset GLO, GAL and BDS by 100, 200 and 300 indexes to keep a single list
% put the two furthest bands in row 1 and 3 since they are preferred
offsetGLO = 100;
offsetGAL = 200;
offsetBDS = 300;

for i = 1:size(obsSolutions,2)
    if isstring(obsSolutions(i).channel)
        channel = obsSolutions(i).channel;
    else
        channel = obsSolutions(i).channel(1:2);
    end

    %--- GPS
    if any(strfind(obsSolutions(i).SV,'GPS')) || contains(obsSolutions(i).SV,'G') 
        if strcmp(channel,'L1') || contains(channel,'G1')
            prnFullIdx(1, obsSolutions(i).PRN) = i;
            obsSolutions(i).freq = 1575.42e6;
            obsSolutions(i).band = 'L1';
        elseif strcmp(channel,'L2') || contains(channel,'G2')
            prnFullIdx(2, obsSolutions(i).PRN) = i;
        elseif strcmp(channel,'L5') || contains(channel,'G5')
            prnFullIdx(3, obsSolutions(i).PRN) = i;
            obsSolutions(i).freq = 1176.45e6;
            obsSolutions(i).band = 'L5';
        end
    end

    %--- GLONASS
    if strfind(obsSolutions(i).SV,'GLO') 
        if strcmp(channel,'G1') || strcmp(channel,'G4')
            prnFullIdx(1, obsSolutions(i).PRN + offsetGLO) = i;
        elseif strcmp(channel,'G2') || strcmp(channel,'G6') 
            prnFullIdx(2, obsSolutions(i).PRN + offsetGLO) = i;
        elseif strcmp(channel,'G3')
            prnFullIdx(3, obsSolutions(i).PRN + offsetGLO) = i;
        end
    end

    %--- Galileo
    if any(strfind(obsSolutions(i).SV,'GAL')) || contains(obsSolutions(i).SV,'E') 
        if strcmp(channel,'E1') || contains(channel,'E1')
            prnFullIdx(1, obsSolutions(i).PRN + offsetGAL) = i;
            obsSolutions(i).freq = 1575.42e6;
            obsSolutions(i).band = 'L1';
        elseif strcmp(channel,'E6')
            prnFullIdx(2, obsSolutions(i).PRN + offsetGAL) = i;
        elseif strcmp(channel,'E5') || strcmp(channel,'E7') || strcmp(channel,'E8') || contains(channel,'E5')
            prnFullIdx(3, obsSolutions(i).PRN + offsetGAL) = i;
            obsSolutions(i).freq = 1176.45e6;
            obsSolutions(i).band = 'L5';
        end
    end

    %--- BeiDou
    if strfind(obsSolutions(i).SV,'BDS')
        if strcmp(channel,'B1') || strcmp(channel,'B2')
            prnFullIdx(1, obsSolutions(i).PRN + offsetBDS) = i;
        elseif strcmp(channel,'B6')
            prnFullIdx(2, obsSolutions(i).PRN + offsetBDS) = i;
        elseif strcmp(channel,'B5') || strcmp(channel,'B7') || strcmp(channel,'B8')
            prnFullIdx(3, obsSolutions(i).PRN + offsetBDS) = i;
        end
    end
end


for i = 1:prnLen
    if sum(~isnan(prnFullIdx(:,i))) >= 2
        %--- Consider the three possible combinations of frequency
        % prefer 1 and 3 if available since they are further
        if ~isnan(prnFullIdx(1,i))
            prnIdx(1,i) = prnFullIdx(1,i);
            if ~isnan(prnFullIdx(3,i))
                prnIdx(2,i) = prnFullIdx(3,i);
            else
                prnIdx(2,i) = prnFullIdx(2,i);
            end
        else
            prnIdx(1,i) = prnFullIdx(2,i);
            prnIdx(2,i) = prnFullIdx(3,i);
        end
    end
end