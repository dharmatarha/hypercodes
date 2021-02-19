function [shiftT, corrV, betas1, Rsq1, est1, betas2, Rsq2, est2, betas3, Rsq3, est3] = couplingPlay(noise)

%% Play around script with the coupling model
%
% We create two loosely correlated random timeseries, shift one of them and
% try to recover that shift using the speaker-listener coupling model of
% Stephens et al, 2010
%
%


%% Parameters, basics

L = 500; % length of vectors we use for this toy example
maxT = 3; % maximum time shift we model
repeatN = 100; % times to re-run the analysis
% noise = 5; % we control the level of correlation by adding noise

% preallocate
shiftT = zeros(repeatN, 1); % time shift in each run (true value)
corrV = zeros(repeatN, 1); % corr values between original time series

% result matrices
betas1 = zeros(2*maxT+1, repeatN);
Rsq1 = zeros(repeatN, 1);
est1 = zeros(repeatN, 1);

betas2 = zeros(2*maxT+1, repeatN);
Rsq2 = zeros(repeatN, 1);
est2 = zeros(repeatN, 1);

betas3 = zeros(2*maxT+1, repeatN);
Rsq3 = zeros(repeatN, 1);
est3 = zeros(repeatN, 1);


%% Loop to do this all many times

for run = 1:repeatN

    %% Create timeseries

    % original, used for modeling the second
    v1 = rand(L,1).*2-1; % uniform between -1 and 1
    % loosely correlated second
    v2 = (v1 + rand(L,1)*noise)/5; 
    corrV(run) = corr(v1,v2); % save out corr value

    % random time lag
    shiftT(run) = randi(2*maxT+1, 1)-(2*maxT+1);
    % shifted coupled timeseries
    if shiftT(run) < 0
        vTarget = [v2(1-shiftT(run):L);zeros(-shiftT(run),1)];
    else
        vTarget = [zeros(shiftT(run),1);v2(1:end-shiftT(run))];
    end

    % shifted versions of the original, each one is a column vector, from -maxT
    % to +maxT
    vShift = zeros(L, 2*maxT+1); % preallocate
    for i = -maxT:1:maxT
        if i < 0
            vShift(:,i+maxT+1) = [v1(1-i:L);zeros(-i,1)];
        else
            vShift(:,i+maxT+1) = [zeros(i,1);v1(1:L-i)];
        end
    end


    %% Get beta values

    % using one understanding of the method in paper: beta = c^-1 * vShifted*vTarget
    betas1(:, run) = cov(vShift)^(-1)*(vShift'*vTarget);
    betas1(:, run) = betas1(:, run)/L;
    SSE = sum((vTarget - (vShift*betas1(:, run))).^2); 
    Rsq1(run) = 1-(SSE / sum(vTarget.^2));  
    [~, b] = max(betas1(:, run));
    est1(run) = b-(maxT+1);
    
    
    % other understanding: beta is given one-by-one
    for z = 1:2*maxT+1
        covZ = cov(vShift(:,z),vTarget)^(-1);
        betas2(z, run) = covZ(1,1)*dot(vShift(:,z),vTarget); 
    end 
    betas2(:, run) = betas2(:, run)/L;
    SSE = sum((vTarget - (vShift*betas2(:, run))).^2); 
    Rsq2(run) = 1-(SSE / sum(vTarget.^2));  
    [~, b] = max(betas2(:, run));
    est2(run) = b-(maxT+1);        
        

    % simple, fast matlab OLS solution
    betas3(:,run) = vShift\vTarget;
    SSE = sum((vTarget - (vShift*betas3(:, run))).^2); 
    Rsq3(run) = 1-(SSE / sum(vTarget.^2));  
    [~, b] = max(betas3(:, run));
    est3(run) = b-(maxT+1);    
    
    if mod(run, 10) == 0
        disp([char(10), 'Done with run ', num2str(run)]);
    end
    
    
end
    

return







