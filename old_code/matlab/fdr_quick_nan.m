function [pQ, RsqFDR, pFDR, lastK] = fdr_quick_nan(Rsq, p)


%% Get k

sortedP = sort(p);
alpha = 0.10;

% a few variables for FDR
m = length(Rsq);
% egamma = 0.577215;

% start with arbitrary k
k = 1000;


tic;
if sortedP(k) < k/m*alpha
    while (sortedP(k) < k/m*alpha) && toc < 20
        lastK = k;
        k = k+1;
    end
elseif sortedP(k) >= k/m*alpha
    while (sortedP(k) >= k/m*alpha) && toc < 20
        lastK = k;
        if k > 1
            k = k-1;
        else
            break;
        end
    end
end

% got the threshold
pQ =sortedP(lastK);

disp(['Last k: ', num2str(lastK)]);
disp(['pQ: ', num2str(pQ)]);


%% Mask of significant columns

mask = p <= pQ;
pFDR = p; 
pFDR(mask) = 1;
RsqFDR = mask.*Rsq;


return
