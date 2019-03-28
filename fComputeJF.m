function JF = fComputeJF(p,yactive,n_ny,eta)
%FCOMPUTEJF Compute Jacobian of F*eta w.r.t. active terms in eta.
%
%	Usage:
%		JF = fComputeJF(p,yactive,n_ny,eta)
%
%	Description:
%		Computes the Jacobian JF of the active terms in F*eta, where eta is
%		a n_ny x N matrix that contains N samples of n_ny known signals
%		(typically nonlinear terms in states and inputs), and where F is a
%		p x n_ny matrix containing the coefficients of these known signals
%		for all p outputs. The indices of the active terms in F are
%		indicated in yactive.
%
%	Output parameters:
%		JF : p*N x n_ny Jacobian matrix with the derivatives of F*eta
%		     w.r.t. the active terms in F
%
%	Input parameters:
%		p : number of outputs
%       yactive : linear indices (see also sub2ind and ind2sub) of the
%                 active elements in transpose(F)
%       n_ny : number of terms in eta
%       eta : n_ny x N matrix with N samples of n_ny known signals
%
%	Example:
%		F = [1 2 3; 4 5 6]; % F matrix with coefficients
%       eta = [1:10:91; 2:10:92; 3:10:93]; % 10 samples of 3 known signals
%       p = size(F,1); % Number of outputs
%       yactive = [1:4 6]; % Fifth element in F not active (i.e. will not be optimized)
%       n_ny = size(eta,1); % Number of signals in eta
%       JF = fComputeJF(p,yactive,n_ny,eta); % Jacobian w.r.t. active elements in F
%       % => JF = [eta.'       zeros(10,2);
%       %          zeros(10,3) eta([1 3],:).'];
%       %       = [derivative of y1 w.r.t. [F(1,1) F(1,2) F(1,3) F(2,1) F(2,3)];
%       %          derivative of y2 w.r.t. [F(1,1) F(1,2) F(1,3) F(2,1) F(2,3)]];
%
%	Versions:
%		1.0 : November 12, 2015
%		1.1 : November 12, 2015
%           Shorter code.
%           Slower on first run, faster on subsequent runs in particular
%           example.
%       1.2 : April 18, 2016
%           Help updated
%
%	Copyright (c) Vrije Universiteit Brussel – dept. ELEC
%   All rights reserved.
%   Software can be used freely for non-commercial applications only.
%   Disclaimer: This software is provided “as is” without any warranty.
%
%	See also sub2ind, ind2sub

%--------------------------------------------------------------------------
% Version 1.1 = Version 1.2
%   Shorter code.
%   Slower on first run, faster on subsequent runs in particular example.
%--------------------------------------------------------------------------
% {

JF = kron(eye(p),eta.'); % Jacobian w.r.t. all elements in F
JF = JF(:,yactive); % Jacobian w.r.t. all active elements in F

%}

%--------------------------------------------------------------------------
% Version 1.0
%--------------------------------------------------------------------------
%{

JF = zeros(p,length(yactive),length(eta)); % p x (#active elements in F) x #samples
temp_ind = 1;
for i = 1:p
    indices = yactive((yactive >= (i-1)*n_ny+1) & (yactive <= i*n_ny))'; % Indices of active terms in output i
    indices = mod(indices-1,n_ny) + 1; % Corresponding column
    JF(i,temp_ind:temp_ind+length(indices)-1,:) = eta(indices,:); % Jacobian of active terms in output i
    temp_ind = temp_ind + length(indices);
end
JF = permute(JF,[3 1 2]); % #samples x p x (#active elements in F)
JF = reshape(JF,p*length(eta),length(yactive)); % p*(#samples) x (#active elements in F)

%}