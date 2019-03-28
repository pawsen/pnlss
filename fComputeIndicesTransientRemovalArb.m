function [indices, NT] = fComputeIndicesTransientRemovalArb(T2,N,p)
%FCOMPUTEINDICESTRANSIENTREMOVALARB Remove transients from arbitrary data.
%
%	Usage:
%		[indices, NT] = fComputeIndicesTransientRemovalArb(T2,N,p)
%
%	Description:
%		Computes the indices to be used with a N x p matrix containing p
%		output signals of length N, such that y(indices) contains the
%		transient-free output(s) of length NT stacked on top of each other
%		(if more than one output). The transient samples to be removed are
%		specified in T2 (T2 = 1:T2 if T2 is scalar).
%
%	Output parameters:
%		indices : vector of indices, such that y(indices) contains the
%	 	          output(s) without transients. If more than one output (p
%	 	          > 1), then y(indices) stacks the transient-free outputs
%	 	          on top of each other.
%       NT : length of the signal without transients
%
%	Input parameters:
%		T2 : scalar indicating how many samples from the start are removed
%		     or vector indicating which samples are removed
%       N : length of the total signal
%       p : number of outputs
%
%	Examples:
%       % One output, T2 scalar
%		N = 1000; % Total number of samples
%       T2 = 200; % First 200 samples should be removed after filtering
%       p = 1; % One output
%       [indices, NT] = fComputeIndicesTransientRemovalArb(T2,N,p);
%       % => indices = (201:1000).'; % Indices of the transient-free output (in uint32 format in version 1.0)
%       % => NT = 800; % Number of samples in the transient-free output
%
%       % Two outputs, T2 scalar
%		N = 1000; % Total number of samples
%       T2 = 200; % First 200 samples should be removed after filtering
%       p = 2; % Two outputs
%       [indices, NT] = fComputeIndicesTransientRemovalArb(T2,N,p);
%       % => indices = ([201:1000 1201:2000]).'; % Indices of the transient-free outputs (in uint32 format in version 1.0)
%       % => NT = 800; % Number of samples in each transient-free output
%       % If y = [y1 y2] is a 1000 x 2 matrix with the two outputs y1 and y2,
%       % then y(indices) = [y1(201:1000);
%       %                    y2(201:1000)]
%       % is a vector with the transient-free outputs stacked on top of
%       % each other
%
%       % One output, T2 is a vector
%       N1 = 1000; % Number of samples in a first data set
%       N2 = 500; % Number of samples in a second data set
%       N = N1 + N2; % Total number of samples
%       T2_1 = 1:200; % Transient samples in first data set
%       T2_2 = 1:100; % Transient samples in second data set
%       T2 = [T2_1 (N1+T2_2)]; % Transient samples
%       p = 1; % One output
%       [indices, NT] = fComputeIndicesTransientRemovalArb(T2,N,p);
%       % => indices = ([201:1000 1101:1500]).'; (in uint32 format in version 1.0)
%       % => NT = 1200;
%
%	Versions:
%		1.0 : November 10, 2015
%       1.1 : April 18, 2016
%           Help updated
%           Removed unsigned 32-bit integer format of indices that was used in the C-code
%
%	Copyright (c) Vrije Universiteit Brussel – dept. ELEC
%   All rights reserved.
%   Software can be used freely for non-commercial applications only.
%   Disclaimer: This software is provided “as is” without any warranty.

%--------------------------------------------------------------------------
% Version 1.1
%   Help updated
%   Removed unsigned 32-bit integer format of indices that was used in the C-code
%--------------------------------------------------------------------------
% {

if length(T2) == 1 % If T2 is scalar
    T2 = (1:T2).'; % Remove all samples up to T2
end

% Which samples are taken into account:
without_T2 = (1:N).'; % All samples in the signal, ...
without_T2(T2) = []; % ... except those that should be removed

NT = length(without_T2); % Length of the transient-free signal(s)
if p > 1 % If more than one output
    indices = zeros(p*NT,1); % Preallocate
    for i = 1:p % Loop over all outputs
        indices((i-1)*NT+1:i*NT) = without_T2 + (i-1)*N; % Stack indices for each output on top of each other
    end
else % Only one output
    indices = without_T2; % Set indices
end

%}

%--------------------------------------------------------------------------
% Version 1.0
%--------------------------------------------------------------------------
%{

if length(T2) == 1 % If T2 is scalar
    T2 = (1:T2).'; % Remove all samples up to T2
end

% Which samples are taken into account:
without_T2 = (1:N).'; % All samples in the signal, ...
without_T2(T2) = []; % ... except those that should be removed

NT = length(without_T2); % Length of the transient-free signal(s)
if p > 1 % If more than one output
    indices = zeros(p*NT,1,'uint32'); % Preallocate + set class of output to unsigned 32-bit integers
    for i = 1:p % Loop over all outputs
        indices((i-1)*NT+1:i*NT) = without_T2 + (i-1)*N; % Stack indices for each output on top of each other
    end
else % Only one output
    indices = uint32(without_T2); % Set class of output to unsigned 32-bit integers
end

%}