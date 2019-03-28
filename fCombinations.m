function out = fCombinations(n,degrees)
%FCOMBINATIONS Lists all nonlinear terms in a multivariate polynomial.
%
%	Usage:
%		out = fCombinations(n,degrees)
%
%	Description:
%		Lists the exponents of all possible monomials in a multivariate
%		polynomial with n inputs. Only the nonlinear degrees in degrees are
%		considered.
%
%	Output parameters:
%		out : ncomb x n matrix of exponents
%
%	Input parameters:
%		n : number of inputs
%       degrees : vector with the degrees of nonlinearity
%
%	Example:
%		% A polynomial with all possible quadratic and cubic terms in the
%		% variables x and y contains the monomials x*x, x*y, y*y, x*x*x,
%		% x*x*y, x*y*y, and y*y*y.
%       out = fCombinations(2,[2 3])
%       % => out = [2 0;        -> x^2 * y^0 = x*x
%       %           1 1;        -> x^1 * y^1 = x*y
%       %           0 2;        -> x^0 * y^2 = y*y
%       %           3 0;        -> x^3 * y^0 = x*x*x
%       %           2 1;        -> x^2 * y^1 = x*x*y
%       %           1 2;        -> x^1 * y^2 = x*y*y
%       %           0 3]        -> x^0 * y^3 = y*y*y
%       % Element (i,j) of out indicates the power to which variable j is
%       % raised in monomial i. For example, out(5,:) = [2 1], which means
%       % that the fifth monomial is equal to x^2*y^1 = x*x*y.
%
%	Versions:
%		1.0 : September 3, 2015
%       1.1 : April 18, 2016
%           Help updated
%
%	Copyright (c) Vrije Universiteit Brussel – dept. ELEC
%   All rights reserved.
%   Software can be used freely for non-commercial applications only.
%   Disclaimer: This software is provided “as is” without any warranty.
%
%	See also nchoosek, lfCombinations

%--------------------------------------------------------------------------
% Version 1.0 = Version 1.1
%--------------------------------------------------------------------------
% {

% Consider only nonlinear degrees
degrees = degrees(degrees > 1);

% Determine total number of combinations
ncomb = 0;
for i = 1:length(degrees)
	ncomb = ncomb + nchoosek(n+degrees(i)-1,degrees(i));
end

% List the exponents of each input in all monomials
out = zeros(ncomb,n); % Preallocating
index = 0; % Running index indicating the last used row in out
for i = 1:length(degrees) % Loop over all nonlinear degrees
	comb_i = lfCombinations(n,degrees(i)); % All combinations in homogeneous polynomial of degree degrees(i)
	out(index + (1:size(comb_i,1)),:) = comb_i; % Insert in out
	index = index + size(comb_i,1); % Update running index
end
end

function combinations = lfCombinations(n,degree)
%LFCOMBINATIONS Lists the exponents of all possible terms in a homogeneous polynomial
    % Monomial representation, e.g. [1 1 2] represents x1*x1*x2
	ncomb = nchoosek(n+degree-1,degree); % Number of combinations in homogeneous polynomial
	monomials = ones(ncomb,degree); % Preallocating, and start from all ones => x1*x1*x1
    for i = 2:ncomb % Loop over all terms
        monomials(i,:) = monomials(i-1,:); % Copy previous monomial
        j = degree; % Index indicating which factor to change
        while monomials(i,j) == n % Try to increase the last factor, ...
            j = j - 1; % ... but if this is not possible, look the previous one that can be increased
        end
        monomials(i,j) = monomials(i,j) + 1; % Increase factor j w.r.t previous monomial, e.g. x1*x1*x1 -> x1*x1*x2
        monomials(i,j+1:degree) = monomials(i,j); % Monomial after x1*x1*xmax is x1*x2*x2, and not x1*x2*xmax
    end

    % Exponents representation, e.g. [2 1] represents x1^2*x2 = x1*x1*x2
	combinations = zeros(ncomb,n); % Preallocating
	for i = 1:size(combinations,1) % Loop over all terms
		for j = 1:n % Loop over all inputs
			combinations(i,j) = sum(monomials(i,:) == j); % Count the number of appearances of input j in monomial i 
		end
	end
end

%}