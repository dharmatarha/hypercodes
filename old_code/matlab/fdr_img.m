function [bFDR, RsqFDR, pFFDR, pQ, lastK] = fdr_img(b, pF, Rsq, alpha)

%% Function for FDR correcting fMRI images
%
% Inputs are from the couplingFMRI function output
%
% Outputs are also in image array format, but values not surviving the FDR
% correction are set to 0 (or to 1 in the case of pFFDR).
% pQ is the highest p-value still considered significant: p(k)
% k is the index we were looking for
%

%% Input checks

if nargin ~= 4
    error([char(10), 'Need four input args, b, pF, Rsq and alpha']);
end
if ~ismatrix(b)
    error([char(10), 'Input arg b is not a matrix']);
end
if ~isvector(pF) || ~isvector(Rsq)
    error([char(10), 'Input arg pF or Rsq is not a vector']);
end
if ~ismembertol(alpha, [0:0.01:1], 0.0001)
    error([char(10), 'Input arg alpha is not between 0 - 1']);
end 


%% Get k

sortedP = sort(pF);

% use Benjamini–Hochberg–Yekutieli procedure for arbitrary dependence
m = length(Rsq);
egamma = 0.577215;

% start with arbitrary k
k = 1000;

% loop until k is found
tic;
if sortedP(k) < k/(m*(log(m)+egamma+(1/2*m)))*alpha
    while (sortedP(k) < k/(m*(log(m)+egamma+(1/2*m)))*alpha) && toc < 20
        lastK = k;
        k = k+1;
    end
elseif sortedP(k) >= k/(m*(log(m)+egamma+(1/2*m)))*alpha
    while (sortedP(k) >= k/(m*(log(m)+egamma+(1/2*m)))*alpha) && toc < 20
        lastK = k;
        if k > 1
            k = k-1;
        else
            break
        end
    end
end

% got the threshold
pQ =sortedP(lastK);


%% Mask of significant columns

mask = pF <= pQ;
pFFDR = pF; 
pFFDR(mask) = 1;
RsqFDR = mask.*Rsq;
bFDR = b; 
for i = 1:m
    if mask(i) == 0
        bFDR(:,i) = zeros(size(b,1),1);
    end
end


return
