function indices = fComputeIndicesTransientRemoval(T1,N,p)
%FCOMPUTEINDICESTRANSIENTREMOVAL Computes indices for transient handling for periodic signals after filtering.
%
%	Usage:
%		indices = fComputeIndicesTransientRemoval(T1,N,p)
%
%	Description:
%       Let u be a vector of length N containing (several realizations of)
%       a periodic signal. Let uTot be a vector containing the signal(s) in
%       u with T1(1) transient points prepended to each realization (see
%       fComputeIndicesTransient). The starting samples of each
%		realization can be specified in T1(2:end). Let yTot be a
%		vector/matrix containing the p outputs of a PNLSS model after
%		applying the input uTot. Then fComputeIndicesTransientRemoval
%		computes the indices to be used with the vectorized form of yTot
%		such that the transient samples are removed from yTot, i.e. y =
%		yTot(indices) contains the steady-state output(s) stacked on top of
%		each other.
%
%	Output parameters:
%		indices : If uTot is a vector containing (several realizations of)
%		          a periodic signal to which T1(1) transient points were
%		          added before each realization, and if yTot is the
%		          corresponding output vector (or matrix if more than one
%		          output), then indices is such that the transient points
%		          are removed from y = yTot(indices).
%                 If p > 1, then indices is a column vector and y =
%                 yTot(indices) is a column vector with the steady state
%                 outputs stacked on top of each other.
%
%	Input parameters:
%		T1 : vector that indicates how the transient is handled. The first
%		     element T1(1) is the number of transient samples that were
%		     prepended to each realization. The other elements T1(2:end)
%		     indicate the starting sample of each realization in the input
%		     signal. If T1 has only one element, T1(2) is put to one.
%       N : length of the input signal containing all realizations
%       p : number of outputs
%
%	Example:
%       Npp = 1000; % Number of points per period
%       R = 2; % Number of phase realizations
%       T = 100; % Number of transient samples
%		T1 = [T 1:Npp:(R-1)*Npp+1]; % Transient handling vector
%       N = R*Npp; % Total number of samples
%       indices_tot = fComputeIndicesTransient(T1,N);
%       % => indices_tot = [901:1000 1:1000 1901:2000 1001:2000]
%       %                = [transient samples realization 1, ...
%       %                   realization 1, ...
%       %                   transient samples realization 2, ...
%       %                   realization 2]
%       p = 1; % One output
%       indices_removal = fComputeIndicesTransientRemoval(T1,N,p);
%       % => indices_removal = [101:1100 1201:2200].'
%       % => indices_tot(indices_removal) = 1:2000
%       %                                 = [realization 1, realization 2]
%
%       p = 2; % More than one output
%       indices_removal = fComputeIndicesTransientRemoval(T1,N,p);
%       % => indices_removal = [101:1100 1201:2200 2301:3300 3401:4400].'
%       % Let u be a vector containing [input realization 1;
%       %                               input realization 2],
%       % then uTot = u(indices_tot) is a vector containing
%       %             [transient samples realization 1;
%       %              input realization 1;
%       %              transient samples realization 2;
%       %              input realization 2]
%       % Let y1 be a vector containing the first output and y2 be a vector
%       % containing the second output when applying uTot as an input to a
%       % PNLSS model, and let yTot = [y1 y2] be a 2200 x 2 matrix with y1
%       % and y2 in its first and second column, respectively.
%       % Note that y1 = yTot(1:2200).' and y2 = yTot(2201:4400).' (see
%       % also ind2sub and sub2ind)
%       % Then yTot(indices_removal) = [y1(101:1100);
%       %                               y1(1201:2200);
%       %                               y2(101:1100);
%       %                               y2(1201:2200)]
%       %                            = [output 1 corresponding to input realization 1;
%       %                               output 1 corresponding to input realization 2;
%       %                               output 2 corresponding to input realization 1;
%       %                               output 2 corresponding to input realization 2]
%
%	Versions:
%		1.0 : November 9, 2015
%       1.1 : April 18, 2016
%           Help updated
%           Removed unsigned 32-bit integer format of indices that was used in the C-code
%
%	Copyright (c) Vrije Universiteit Brussel – dept. ELEC
%   All rights reserved.
%   Software can be used freely for non-commercial applications only.
%   Disclaimer: This software is provided “as is” without any warranty.
%
%	See also fComputeIndicesTransient, ind2sub, sub2ind, fVec

%--------------------------------------------------------------------------
% Version 1.1
%   Help updated
%   Removed unsigned 32-bit integer format of indices that was used in the C-code
%--------------------------------------------------------------------------
% {

if length(T1) == 1% If starting samples of realizations not specified, ...
	T1 = [T1 1]; % ... then one realization that starts at first sample
end

Trans = T1(1); % Number of transient points
T = [fVec(T1(2:end)); N+1]'; % Vector containing all starting samples of the realizations and the length of the total signal plus one
indices = []; % Initialize indices as an empty matrix
for i = 1:length(T)-1 % Loop over all realizations
    indices = [indices (T(i):T(i+1)-1)+i*Trans]; % Concatenate indices without transient samples
end

if p > 1 % If more than one output
    NT = N + T1(1)*(length(T1)-1); % Total number of samples per output = number of samples without + with transients
    temp = zeros(p*N,1); % Preallocate
    for i = 1:p % Loop over all outputs
        temp((i-1)*N+1:i*N) = indices + (i-1)*NT; % Stack indices without transient samples on top of each other
    end
    indices = temp; % Collect result from temporary variable
end

%}

%--------------------------------------------------------------------------
% Version 1.0
%--------------------------------------------------------------------------
%{

if length(T1) == 1% If starting samples of realizations not specified, ...
	T1 = [T1 1]; % ... then one realization that starts at first sample
end
   

Trans = T1(1); % Number of transient points
T = [fVec(T1(2:end)); N+1]'; % Vector containing all starting samples of the realizations and the length of the total signal plus one
indices = uint32([]); % Set class of output to unsigned 32-bit integers
for i = 1:length(T)-1 % Loop over all realizations
    indices = [indices (T(i):T(i+1)-1)+i*Trans]; % Concatenate indices without transient samples
end

if p > 1 % If more than one output
    NT = N + T1(1)*(length(T1)-1); % Total number of samples per output = number of samples without + with transients
    temp = zeros(p*N,1,'uint32'); % Preallocate + set class of output to unsigned 32-bit integers
    for i = 1:p % Loop over all outputs
        temp((i-1)*N+1:i*N) = indices + (i-1)*NT; % Stack indices without transient samples on top of each other
    end
    indices = temp; % Collect result from temporary variable
end

%}