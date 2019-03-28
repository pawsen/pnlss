function out = fJNL(input,Edwdx,C,Fdwdx,active)
%FJNL Compute Jacobian w.r.t. A, B, and E by filtering an alternative state-space model.
%
%	Usage:
%		out = fJNL(input,Edwdx,C,Fdwdx,active)
%
%	Description:
%		Computing the Jacobian of the output y(t) of a nonlinear
%		state-space model
%
%           x(t+1) = A x(t) + B u(t) + E zeta(x(t),u(t))
%           y(t)   = C x(t) + D u(t) + F eta(x(t),u(t))
%
%       w.r.t. the elements in the A, B, and E matrices can be performed by
%       filtering an alternative nonlinear state-space model.
%       Let JA(i,j)(t) be the partial derivative of y(t) w.r.t. A(i,j) and
%       let xA(i,j)(t) be the partial derivative of x(t) w.r.t. A(i,j).
%       Similarly, consider JB(i,j)(t), JE(i,j)(t), xB(i,j)(t), and
%       xE(i,j)(t) to be the partial derivative of y(t) and x(t) w.r.t.
%       B(i,j) and E(i,j).
%       Then the Jacobians JA, JB, and JE can be computed by filtering
%       three state-space models (Paduart, 2008):
%
%           xA(i,j)(t+1) = A xA(i,j)(t) + I(i,j) x(t) + E partial(zeta(t))/partial(x(t)) xA(i,j)(t) 
%           JA(i,j)(t)   = C xA(i,j)(t)               + F partial(eta(t))/partial(x(t))  xA(i,j)(t) 
%
%           xB(i,j)(t+1) = A xB(i,j)(t) + I(i,j) u(t) + E partial(zeta(t))/partial(x(t)) xB(i,j)(t) 
%           JB(i,j)(t)   = C xB(i,j)(t)               + F partial(eta(t))/partial(x(t))  xB(i,j)(t) 
%
%           xE(i,j)(t+1) = A xE(i,j)(t) + I(i,j) zeta(t) + E partial(zeta(t))/partial(x(t)) xE(i,j)(t) 
%           JE(i,j)(t)   = C xE(i,j)(t)                  + F partial(eta(t))/partial(x(t))  xE(i,j)(t) 
%
%       where I(i,j) is a matrix of the appropriate dimensions with all
%       zeros, except in position (i,j), where I(i,j) = 1 (see also fOne).
%       fJNL computes the Jacobian w.r.t. either A, B, or E, depending on
%       the provided input (x(t), u(t), or zeta(t)).
%
%	Output parameters:
%		out : p x N x nactive Jacobian matrix with the partial derivatives
%		      of y(t) w.r.t. the active elements in A, B, or E (depending
%		      on the provided input), where p is the number of outputs, N
%		      is the number of samples, and nactive is the number of active
%		      elements in A, B, or E.
%
%	Input parameters:
%		input : N x npar matrix with the input samples in the alternative
%		        state-space model (input = states x(t) and npar = n if JA
%		        is calculated, input = input u(t) and npar = m if JB is
%		        calculated, and input = zeta(t) and npar = n_nx if JE is
%		        calculated).
%		Edwdx : n x n x N matrix with N samples of A + E partial(zeta(t))/partial(x(t))
%		C : p x n matrix C
%		Fdwdx : p x n x N matrix with N samples of F partial(eta(t))/partial(x(t))
%		active : nactive x 1 vector with the linear indices (see also
%		         sub2ind and ind2sub) in A.', B.', or E.' indicating the
%		         active elements of A, B, or E (active elements = elements
%		         on which optimization will be done, see also
%		         fSelectActive).
%
%   Reference:
%       Paduart J. (2008). Identification of nonlinear systems using
%       polynomial nonlinear state space models. PhD thesis, Vrije
%       Universiteit Brussel.
%
%	Versions:
%		1.0 : November 12, 2015
%       1.1 : April 20, 2016
%           Help updated
%
%	Copyright (c) Vrije Universiteit Brussel – dept. ELEC
%   All rights reserved.
%   Software can be used freely for non-commercial applications only.
%   Disclaimer: This software is provided “as is” without any warranty.
%
%	See also fOne, sub2ind, ind2sub, fSelectActive

%--------------------------------------------------------------------------
% Version 1.0 = Version 1.1
%--------------------------------------------------------------------------
% {

[p,n] = size(C); % Number of outputs and number of states
[N,npar] = size(input); % Number of samples and number of inputs in alternative state-space model
nactive = length(active); % Number of active parameters in A, B, or E

Fdwdx = Fdwdx + repmat(C,[1,1,N]); % Add C to F partial(eta(t))/partial(x(t)) for all samples at once (faster than doing this inside the loop)

out = zeros(p,N,nactive); % Preallocate
for k = 1:nactive % Loop over all active elements in A, B, or E
    j = mod(active(k)-1,npar) + 1; % Which column in A, B, or E matrix
    i = (active(k)-j)/npar + 1; % Which row in A, B, or E matrix
    Jprev = zeros(n,1); % Initialize partial derivative of x(0) w.r.t. A(i,j), B(i,j), or E(i,j)
    for t = 2:N % Loop over all samples
        % Calculate state update alternative state-space model at time t 
        J = Edwdx(:,:,t-1)*Jprev; % Terms in alternative states at time t-1
        J(i) = J(i) + input(t-1,j); % Term in alternative input at time t-1
        % Calculate output alternative state-space model at time t
        out(:,t,k) = Fdwdx(:,:,t)*J;
        % Update previous state alternative state-space model
        Jprev = J;
    end
end

%}

%about as slow on a full model (i.e. all terms active)
%{

[p,n] = size(C); % Number of outputs and number of states
[N,npar] = size(input); % Number of samples and number of inputs in alternative state-space model
nactive = length(active); % Number of active parameters in A, B, or E

Fdwdx = Fdwdx + repmat(C,[1,1,N]); % Add C to F partial(eta(t))/partial(x(t)) for all samples at once (faster than doing this inside the loop)

i_active = 1; % Index of active terms in A, B, or E
out = zeros(p,N,nactive); % Preallocate
for i = 1:n % Loop over all states
    for j = 1:npar % Loop over all inputs in alternative state-space model
        if active(i_active) == (i-1)*npar + j % If A(i,j), B(i,j), or E(i,j) is an active parameter
            Jprev = zeros(n,1); % Initialize partial derivative of x(0) w.r.t. A(i,j), B(i,j), or E(i,j)
            for t = 2:N % Loop over all samples
                % Calculate state update alternative state-space model at time t 
                J = Edwdx(:,:,t-1)*Jprev; % Terms in alternative states at time t-1
                J(i) = J(i) + input(t-1,j); % Term in alternative input at time t-1
                % Calculate output alternative state-space model at time t
                out(:,t,i_active) = Fdwdx(:,:,t)*J;
                % Update previous state alternative state-space model
                Jprev = J;
            end
            i_active = i_active + 1; % Update index of active terms
        end
    end
end

%}