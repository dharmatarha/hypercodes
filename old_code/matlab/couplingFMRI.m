function [b, Rsq, F, pF, pP] = couplingFMRI(speaker, listener, maxT, stats)

%% Function to calculate ISC between speaker-listener
%
% [b, Rsq, F, pF, pP] = couplingFMRI(speaker, listener, maxT, stats=0)
%
% Based on Stephens et al, 2010 and Silbert et al, 2014
%
% Input args
% speaker:      Model data, epi timeseries as matrix (TR x voxels)
% listener:     Data to model, epi timeseries as matrix (TR x voxels)
% maxT:         Maximum lag to include in model, given as no. of TR
% stats:        If set, the function also performs a random permutation
%               test for each voxel. The number of permutations is given by
%               the input argument, e.g., if stats == 1000, 1000
%               permutations will be done for each voxel timeseries.
%               Default is no permutation test.
%
% Outputs:
% b:            Beta vectors, weights for all shifted timeseries, for all
%               voxels, a matrix sized (2*maxT+1 x voxels) 
% Rsq:          Overall fit for each voxel, Rsquared, in row vector
% F:            F statistic for the model at each voxel, row vector
% pF:           P-values for the F statistic at each voxel, row vector 
% pP:           P-values for each voxel, as given by the permutation test
%
% Notes:
% - Change random permutation to phase-scrambling permutation, preserving
% the power spectrum of the original time-series at each permutation
% - For permutation test we take Rsquared as our test statistic - should we
% choose stg else?
% - Optimize?
% - Double-check the F-statistic part, especially the DFs


%% Input checks

if nargin < 3
    error([char(10), 'Input args speaker, listener and maxT are needed']);
end
if nargin == 3
    stats = 0;
end
if ~ismatrix(speaker) || ~ismatrix(listener)
    error([char(10), 'Input arg speaker or listener is not a matrix']);
end
if size(speaker, 1) > size(speaker, 2) || size(listener, 1) > size(listener, 2)
    error([char(10), 'Input arg speaker or listener has wrong shape, more rows than columns']);
end
if ~isequal(size(speaker), size(listener))
    error([char(10), 'Input args speaker and listener should have same size']);
end
if ~ismember(maxT, 1:10)
    error([char(10), 'Input arg maxT should be int between 1 - 10']);
end
if (stats < 0) || (stats > 10000)
    error([char(10), 'Input arg stats should be between 0 - 10000, eh?']);
end

disp([char(10), 'Started couplingFMRI.m with inputs ', inputname(1), ...
    ' (speaker epi matrix), ', inputname(2), ' (listener epi matrix), ', char(10), ...
    num2str(maxT), ' (max time lag) and ', num2str(stats), ' (stats param)']);


%% Basics

% no. of voxels
voxelN = size(speaker, 2);
% no. of TRs
trN = size(speaker, 1);
% preallocate results matrices
b = zeros(2*maxT+1, voxelN);
Rsq = zeros(1, voxelN);
F = zeros(1, voxelN);
pF = zeros(1, voxelN);
pP = zeros(1, voxelN);

disp([char(10), 'Data dimensions: ', num2str(trN), ' X ', num2str(voxelN)]);
disp('Modelling listener data using lagged time series from speaker...');


%% Loop through voxels

for vox = 1:voxelN
    
    % create model matrix from time-shifted vectors of speaker
    vShift = zeros(trN, 2*maxT+1); % preallocate
    for i = -maxT:1:maxT
        if i < 0
            vShift(:,i+maxT+1) = [speaker(1-i:trN,vox);zeros(-i,1)];
        else
            vShift(:,i+maxT+1) = [zeros(i,1);speaker(1:trN-i,vox)];
        end
    end

    % solve linear model vShift*b=y where y is the vector from listener's
    % data
    b(:,vox) = vShift\listener(:,vox);
    
    % get Rsquared
    SSE = sum((listener(:,vox)-(vShift*b(:,vox))).^2); 
    Rsq(vox) = 1-(SSE / sum(listener(:,vox).^2));
    % get F statistic for model
    % SSM = sum((vShift*b).^2);
    SSM = sum(listener(:,vox).^2) - SSE;
    DFM = 2*maxT; % +1 ?
    DFE = trN-(2*maxT+1); % -1 ?
    F(vox) = (SSM/DFM)/(SSE/DFE);
    % get p-value for F stat
    pF(vox) = 1-fcdf(F(vox), DFM, DFE);
    
    
    %% Permutation part in voxel loop
    if stats
        
        % preallocate permutation results 
        permRsq = zeros(stats,1);
        
        % prepare fft of voxel time series in advance
        sigY = fft(speaker(:, vox));
        [sigP, sigMagn] = cart2pol(real(sigY), imag(sigY)); % get fft in polar coordinates
        phases = sigP(2:ceil(trN/2));
        
        for perm = 1:stats

            % phase-randomized permutation of speaker voxel time series
            pRandom = phases(randperm(length(phases))); % randomize phase component of fft
            % keep conjugate symmetry so that we have real data at the end
            if mod(trN, 2) == 0
                pNew = [sigP(1); pRandom; sigP(trN/2+1); flipud(pRandom.*(-1))];
            else
                pNew = [sigP(1); pRandom; flipud(pRandom.*(-1))];
            end
            % get complex from polar coordinates
            [xCoeff, yCoeff] = pol2cart(pNew, sigMagn);
            sigYNew = xCoeff + yCoeff*i;
            % transform back into time domain
            permV = ifft(sigYNew, 'symmetric'); % 'symmetric' ensures that rounding errors do not result in complex time series
            
            % sanity check
            if ~isreal(permV)
                error('Phase randomization failed during permutation, permuted signal was not real array, but probably complex');
            end
            
            % create shifted time series matrix
            permVShift = zeros(trN, 2*maxT+1); % preallocate
            for i = -maxT:1:maxT
                if i < 0
                    permVShift(:,i+maxT+1) = [permV(1-i:trN);zeros(-i,1)];
                else
                    permVShift(:,i+maxT+1) = [zeros(i,1);permV(1:trN-i)];
                end
            end

            % solve linear model permVShift*b=y
            permB = vShift\listener(:,vox);        
            % get Rsquared
            permSSE = sum((listener(:,vox)-(permVShift*permB)).^2); 
            permRsq(perm) = 1-(permSSE / sum(listener(:,vox).^2));        
        end
        
        % get p-value from permutated Rsq values
        pP(vox) = sum(permRsq >= Rsq(vox))/stats;
        
    end
        
    % occassional feedback
    if mod(vox, 10000) == 0
        disp([char(10), 'Finished with voxel no. ', num2str(vox)]);
    end
    
    
end


return










