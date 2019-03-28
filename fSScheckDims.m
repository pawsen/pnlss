function [n,m,p] = fSScheckDims(A,B,C,D)
%FSSCHECKDIMS Check consistency of state-space dimensions.
%
%	Usage:
%		[n,m,p] = fSScheckDims(A,B,C,D)
%
%	Description:
%		Returns the number of states, inputs, and outputs of a linear
%		state-space model and checks if the sizes of the state-space
%		matrices are consistent. Produces an error message if they are not
%		consistent.
%
%	Output parameters:
%		n : model order
%       m : number of inputs
%       p : number of outputs
%
%	Input parameters:
%		A : n x n state matrix
%       B : n x m input matrix
%       C : p x n output matrix
%       D : p x m feed-through matrix
%
%	Versions:
%		1.0 : September 4, 2015
%       1.1 : April 20, 2016
%           Help updated
%       1.2 : August 12, 2016
%           Help updated
%
%	Copyright (c) Vrije Universiteit Brussel – dept. ELEC
%   All rights reserved.
%   Software can be used freely for non-commercial applications only.
%   Disclaimer: This software is provided “as is” without any warranty.

%--------------------------------------------------------------------------
% Version 1.0 = Version 1.1 = Version 1.2
%--------------------------------------------------------------------------
% {

% Determine sizes of A, B, C, and D matrix
[n1,n2] = size(A);
[n3,m1] = size(B);
[p1,n4] = size(C);
[p2,m2] = size(D);

% Check whether the number of states is consistent
if (n1 == n2) && (n2 == n3) && (n3 == n4)
    n = n1;
else
    error('Inconsistent dimensions');
end

% Check whether the number of inputs is consistent
if m1 == m2
    m = m1;
else
    error('Inconsistent dimensions');
end

% Check whether the number of outputs is consistent
if p1 == p2
    p = p1;
else
    error('Inconsistent dimensions');
end

%}