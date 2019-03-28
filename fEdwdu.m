function out = fEdwdu(contrib,pow,coeff,E,nx,m)
%FEDWDU Multiply a matrix E with the derivative w.r.t. u of a polynomial w(x,u).
%
%	Usage:
%		out = fEdwdu(contrib,pow,coeff,E,nx,n)
%
%	Description:
%		Multiplies a matrix E with the derivative of a polynomial w(x,u)
%		w.r.t. the m elements in u. The samples of x and u are in a vector
%		contrib. The derivative of w(x,u) w.r.t. u is given by the
%		exponents in x and u (given in pow) and the corresponding
%		coefficients (given in coeff). The maximum degree of a variable (an
%		x or a u) in w(x,u) is given in nx.
%
%	Output parameters:
%		out : n_out x m x N matrix that is the product of E and the
%		      derivative of the polynomial w(x,u) w.r.t. the elements in u
%		      at all samples.
%
%	Input parameters:
%		contrib : (n+m) x N matrix containing N samples of the signals x
%		          and u
%       pow : n_nx x (n+m) x (n+m) containing the exponents of the
%             derivatives of w(x,u) w.r.t. x and u, i.e. pow(i,j,k)
%             contains the exponent of contrib j in the derivative of the
%             ith monomial w.r.t. contrib k.
%       coeff : n_nx x (n+m) matrix containing the corresponding
%               coefficients, i.e. coeff(i,k) contains the coefficient of
%               the derivative of the ith monomial in w(x,u) w.r.t. contrib
%               k.
%       E : n_out x n_nx matrix
%       nx : maximum degree of a variable (an x or a u) in w(x,u)
%       m : number of u signals w.r.t. which derivatives are taken
%
%	Example:
%       % Consider w(x1,x2,u) = [x1^2    and E = [1 3 5
%       %                        x1*x2;           2 4 6]
%       %                        x2*u^2]
%       % then the derivatives of E*w w.r.t. x1 and x2 are given by
%       % E*[2*x1 0
%       %    1*x2 1*x1
%       %    0    1*u^2]
%       % and the derivative of w w.r.t. u is given by [0
%       %                                               0
%       %                                               2*x2*u]
%		E = [1 3 5; 2 4 6];
%       pow = zeros(3,3,3);
%       pow(:,:,1) = [1 0 0;
%                     0 1 0;
%                     0 0 0]; % Derivative w.r.t. x1 has terms 2*x1, 1*x2, and 0
%       pow(:,:,2) = [0 0 0;
%                     1 0 0;
%                     0 0 2]; % Derivative w.r.t. x2 has terms 0, 1*x1, and 1*u^2
%       pow(:,:,3) = [0 0 0;
%                     0 0 0;
%                     0 1 1]; % Derivative w.r.t. u has terms 0, 0, and 2*x2*u
%       coeff = [2 0 0;
%                1 1 0;
%                0 1 2];
%       nx = 2; % Maximum second degree factor in monomials of w (x1^2 in first monomial, u^2 in third monomial)
%       m = 1; % One signal u
%       contrib = randn(3,10); % Ten random samples of signals x1, x2, and u
%       out = fEdwdu(contrib,pow,coeff,E,nx,m);
%       % => out(:,:,t) = E*[0
%       %                    0
%       %                    2*contrib(2,t)*contrib(3,t)]
%
%	Versions:
%		1.0 : January 14, 2016
%       1.1 : April 19, 2016
%           Help updated
%
%	Copyright (c) Vrije Universiteit Brussel – dept. ELEC
%   All rights reserved.
%   Software can be used freely for non-commercial applications only.
%   Disclaimer: This software is provided “as is” without any warranty.
%
%   See also fEdwdx

%--------------------------------------------------------------------------
% Version 1.0 = Version 1.1
%--------------------------------------------------------------------------
% {

[n_all,N] = size(contrib); % n_all = number of signals x and u; N = number of samples
n = n_all - m; % Number of signals x
[n_out,n_nx] = size(E); % n_out = number of rows in E; n_nx = number of monomials in w
out = zeros(n_out,m,N); % Preallocate
for k = 1:m % Loop over all signals u w.r.t. which derivatives are taken
    A = coeff(:,n+k)*ones(1,N); % Repeat coefficients of derivative of w w.r.t. u_k
    for j = 1:n_all % Loop over all signals x and u
        for i = 1:n_nx % Loop over all monomials
            A(i,:) = A(i,:).*contrib(j,:).^pow(i,j,n+k); % Derivative of monomial i w.r.t. u_k
        end
    end
    out(:,k,:) = E*A; % E times derivative of w w.r.t. u_k
end

%}