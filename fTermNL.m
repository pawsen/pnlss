function out = fTermNL(contrib,pow,max_degree)
%FTERMNL Construct polynomial terms.
%
%	Usage:
%		out = fTermNL(contrib,pow,max_degree)
%
%	Description:
%		out = fTermNL(contrib,pow,max_degree) computes polynomial terms,
%		where contrib contains the input signals to the polynomial and pow
%		contains the exponents of each term in each of the inputs. The
%		maximum degree of an individual input is given in max_degree.
%
%	Output parameters:
%		out : nterms x N matrix with N samples of each term
%
%	Input parameters:
%		contrib : (n+m) x N matrix with N samples of the input signals to
%		          the polynomial. Typically, these are the n states and the
%		          m inputs of the nonlinear state-space model.
%       pow : nterms x (n+m) matrix with the exponents of each term in
%             each of the inputs to the polynomial.
%       max_degree : maximum degree in an individual input of the
%                    polynomial
%
%	Example:
%		n = 2; % Number of states
%       m = 1; % Number of inputs
%       N = 1000; % Number of samples
%       x = randn(n,N); % States
%       u = randn(m,N); % Input
%       contrib = [x; u]; % States and input combined
%       pow = [2 0 0;
%              1 1 0;
%              1 0 1;
%              0 2 0;
%              0 1 1;
%              0 0 2]; % All possible quadratic terms in states and input: x1^2, x1*x2, x1*u, x2^2, x2*u, u^2 
%       max_degree = max(max(pow)); % Maximum degree in an individual state or input
%       out = fTermNL(contrib,pow,max_degree);
%       % => out = [x(1,:).^2;
%       %           x(1,:).*x(2,:);
%       %           x(1,:).*u;
%       %           x(2,:).^2;
%       %           x(2,:).*u;
%       %           u.^2];
%
%	Versions:
%		1.0 : November 18, 2015
%       1.1 : April 20, 2016
%           Help updated
%
%	Copyright (c) Vrije Universiteit Brussel – dept. ELEC
%   All rights reserved.
%   Software can be used freely for non-commercial applications only.
%   Disclaimer: This software is provided “as is” without any warranty.
%
%   See also fCombinations

%--------------------------------------------------------------------------
% Version 1.0 = Version 1.1
%--------------------------------------------------------------------------
% {

N = size(contrib,2); % Number of samples
nterms = size(pow,1); % Number of terms
out = zeros(nterms,N); % Preallocate
pow = repmat(permute(pow,[2 3 1]),[1 N 1]); % (n+m) x N x nterms
for i = 1:nterms % Loop over all terms
	out(i,:) = prod(contrib.^pow(:,:,i)); % All samples of term i
end

%}