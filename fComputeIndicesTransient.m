function indices = fComputeIndicesTransient(T1,N)
%FCOMPUTEINDICESTRANSIENT Computes indices for transient handling for periodic signals before filtering.
%
%	Usage:
%		indices = fComputeIndicesTransient(T1,N)
%
%	Description:
%		fComputeIndicesTransient computes the indices to be used with a
%		vector u of length N that contains (several realizations of) a
%		periodic signal, such that u(indices) has T1(1) transient samples
%		prepended to each realization. The starting samples of each
%		realization can be specified in T1(2:end).
%       Like this, steady-state data can be obtained from a PNLSS model by
%       using u(indices) as an input signal to a PNLSS model (see
%       fFilterNLSS) and removing the transient samples afterwards (see
%       fComputeIndicesTransientRemoval).
%
%	Output parameters:
%		indices : indices of a vector u that contains (several realizations
%		          of) a periodic signal, such that u(indices) has a number
%		          of transient samples added before each realization
%
%	Input parameters:
%		T1 : vector that indicates how the transient is handled. The first
%		     element T1(1) is the number of transient samples that should
%		     be prepended to each realization. The other elements T1(2:end)
%		     indicate the starting sample of each realization in the
%		     signal. If T1 has only one element, T1(2) is put to one.
%       N : length of the signal containing all realizations
%
%	Example:
%       Npp = 1000; % Number of points per period
%       R = 2; % Number of phase realizations
%       T = 100; % Number of transient samples
%		T1 = [T 1:Npp:(R-1)*Npp+1]; % Transient handling vector
%       N = R*Npp; % Total number of samples
%       indices = fComputeIndicesTransient(T1,N);
%       % => indices = [901:1000 1:1000 1901:2000 1001:2000]
%       %            = [transient samples realization 1, ...
%       %               realization 1, ...
%       %               transient samples realization 2, ...
%       %               realization 2]
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
%	See also fFilterNLSS, fComputeIndicesTransientRemoval, fVec

%--------------------------------------------------------------------------
% Version 1.1
%   Help updated
%   Removed unsigned 32-bit integer format of indices that was used in the C-code
%--------------------------------------------------------------------------
% {

if T1(1) ~= 0
    Trans = T1(1); % Number of transient points
    if length(T1) == 1 % If starting samples of realizations not specified, ...
        T1 = [T1 1]; % ... then one realization that starts at first sample
    end
    T1 = [fVec(T1(2:end)); N+1].'; % Vector containing all starting samples of the realizations and the length of the total signal plus one
    indices = []; % Initialize indices as an empty matrix
    for i = 1:length(T1)-1 % Loop over all realizations
        transientpart = T1(i+1)-1-mod(Trans-1:-1:0,T1(i+1)-T1(i)); % Indices transient part realization i (assuming periodic signal)
        normalpart    = T1(i):T1(i+1)-1; % Indices realization i
        indices = [indices transientpart normalpart]; % Concatenate indices transient and normal part and put together with indices of previous realizations
    end
else
    indices = 1:N; % No transient points => output = all indices of the signal
end

%}

%--------------------------------------------------------------------------
% Version 1.0
%--------------------------------------------------------------------------
%{

if T1(1) ~= 0
    Trans = T1(1); % Number of transient points
    if length(T1) == 1 % If starting samples of realizations not specified, ...
        T1 = [T1 1]; % ... then one realization that starts at first sample
    end
    T1 = [fVec(T1(2:end)); N+1].'; % Vector containing all starting samples of the realizations and the length of the total signal plus one
    indices = uint32([]); % Set class of output to unsigned 32-bit integers
    for i = 1:length(T1)-1 % Loop over all realizations
        transientpart = T1(i+1)-1-mod(Trans-1:-1:0,T1(i+1)-T1(i)); % Indices transient part realization i (assuming periodic signal)
        normalpart    = T1(i):T1(i+1)-1; % Indices realization i
        indices = [indices transientpart normalpart]; % Concatenate indices transient and normal part
    end
else
    indices=uint32(1:N); % No transient points => output = all indices of the signal in unsigned 32-bit integer format
end

%}