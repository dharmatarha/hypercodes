function [n, s, sP, p, magn, pNew, xCoeff, yCoeff, yNew, sNew] = phaseRandom

% phase randomization play

% signal
n = 500; % length
s = rand(n,1).*2-1;

% fft
y = fft(s);

% sanity check
sP = ifft(y);
diff = sum(abs(s-sP));
disp([char(10), 'Back-and-forth Fourier difference was ', num2str(diff)]);

% get phases and magnitude
[p, magn] = cart2pol(real(y), imag(y));

% sanity check for symmetry
if ~isequal(p(2:ceil(n/2)),flipud(p(floor(n/2)+2:end).*(-1)))
    error('Oops, phases have a problem...');
end

% randomize phases
phases = p(2:ceil(n/2));
pRandom = phases(randperm(length(phases)));
if mod(length(s), 2) == 0
    pNew = [p(1); pRandom; p(n/2+1); flipud(pRandom.*(-1))];
else
    pNew = [p(1); pRandom; flipud(pRandom.*(-1))];
end

% get complex from polar coordinates
[xCoeff, yCoeff] = pol2cart(pNew, magn);
yNew = xCoeff + yCoeff*i;

% ifft 
sNew = ifft(yNew, 'symmetric');

