function [Ls, P, meanB] = contrastT(bAll, voxelMask, contrast)

%% Quick and dirty contrast analysis function
%
% 
%

% we need a column vector as contrasts
if size(contrast,2) > size(contrast,1)
    contrast=contrast';
    disp('contrast vector was not a column vector, switched it');
end

voxelN = length(voxelMask);
N = size(bAll,1);
validV = sum(voxelMask);
Ts = size(bAll,2);
Ls =nan(N, voxelN);
P = nan(1, voxelN);
meanB = nan(Ts, voxelN);

for v = 1:voxelN
    
    if voxelMask(v)
        
        meanB (:, v) = squeeze(mean(bAll(:,:,v),1));
        
        Ls(:, v) = squeeze(bAll(:,:,v))*contrast;
        
        [~, P(v)] = ttest(Ls(:, v));
        
    end
    
end

return