function covY = fCovarY(y)
%FCOVARY Compute covariance matrix output spectra due to output noise from time domain output signals.
%
%	Usage:
%		covY = fCovarY(y)
%
%	Description:
%		Calculate frequency domain covariance matrix of the output spectrum
%		starting from time domain output signals.
%
%	Output parameters:
%		covY : p x p x NFD covariance matrix of the output(s). Calculates
%		       variations along the periods and averages over the
%		       realizations.
%
%	Input parameters:
%		y : N x p x R x P output signal(s), where N is the number of
%		    samples, p is the number of outputs, R is the number of
%		    realizations, and P is the number of periods
%
%   Example:
%       y0 = randn(1000,1); % Noise-free signal
%       y0 = repmat(y0,3,1); % Three periods
%       y0 = reshape(y0,[1000,1,1,3]); % N x p x R x P
%       [b,a] = cheby1(2,10,2*0.1); % Coefficients noise filter
%       e = randn(size(y0)); % White noise
%       v = filter(b,a,e); % Colored noise
%       y = y0 + v; % Noisy signal
%       covY = fCovarY(y); % Noise covariance
%       figure
%       plot(db(freqz(b,a,2*pi*(1:500)/1000)*sqrt(1000)/sqrt(3)),'b')
%       hold on
%       plot(db(squeeze(sqrt(covY(:,:,1:500)))),'r')
%       xlabel('Frequency line')
%       ylabel('Amplitude (dB)')
%       title('Noise variance on the averaged signal (average over 3 periods)')
%       legend('True','Estimated')
%
%	Versions:
%		1.0 : November 10, 2015
%       1.1 : April 18, 2016
%           Help updated
%
%	Copyright (c) Vrije Universiteit Brussel – dept. ELEC
%   All rights reserved.
%   Software can be used freely for non-commercial applications only.
%   Disclaimer: This software is provided “as is” without any warranty.
%
%	See also 

%--------------------------------------------------------------------------
% Version 1.0 = Version 1.1
%--------------------------------------------------------------------------
% {

[N,p,R,P] = size(y); % Number of samples, outputs, realizations, and periods
Y = fft(y); % Go to the frequency domain
NFD = floor(N/2); % Number of bins in positive half of the spectrum
Y = Y(1:NFD,:,:,:); % Take only the positive half of the spectum into account (negative half = conjugate of positive half)
Ymean = mean(Y,4); % Average over the periods
NY = Y - repmat(Ymean,[1 1 1 P]); % Variations over periods
NY = permute(NY,[2 3 4 1]); % p x R x P x NFD
covY = zeros(p,p,R,NFD); % Preallocate
for i = 1:NFD % Loop over frequency bins
    for r = 1:R % Loop over realizations
        temp = zeros(p,p); % Initialize temporary sum to zero matrix
        for pp = 1:P % Loop over periods
            temp = temp + fVec(NY(:,r,pp,i))*(fVec(NY(:,r,pp,i)))'; % Sum outer products of variations over periods and its complex conjugate
        end
        covY(:,:,r,i) = temp/(P-1)/P; % Output covariance of realization r at frequency bin i, averaged over P periods
    end
end
covY = mean(covY,3); % Compute average over all realizations
covY = permute(covY,[1 2 4 3]); % p x p x NFD

%}