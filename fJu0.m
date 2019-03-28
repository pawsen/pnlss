function out = fJu0(B_EdwxIdu0,Edwdx,C,Fdwdx,active)
%FJU0 Compute Jacobian w.r.t. u0 by filtering an alternative state-space model.
%
%	Usage:
%		out = fJu0(B_EdwxIdu0,Edwdx,C,Fdwdx,active)
%
%	Description:
%		Computing the Jacobian of the output y(t) of a nonlinear
%		state-space model
%
%           x(t+1) = A x(t) + B u(t) + E zeta(x(t),u(t))
%           y(t)   = C x(t) + D u(t) + F eta(x(t),u(t))
%
%       w.r.t. the elements in the initial input u0 can be performed by
%       filtering an alternative nonlinear state-space model.
%       Let Ju0(i)(t) be the partial derivative of y(t) w.r.t. u0(i) and
%       let xu0(i)(t) be the partial derivative of x(t) w.r.t. u0(i).
%       Then the Jacobian Ju0 can be computed by filtering an alternative
%       state-space model
%
%           xu0(i)(t+1) = A xu0(i)(t) + (B I(i,1) + E partial(zeta(0))/partial(u0(i))) delta(t,0) + E partial(zeta(t))/partial(x(t)) xu0(i)(t) 
%           Ju0(i)(t)   = C xu0(i)(t)                                                             + F partial(eta(t))/partial(x(t))  xu0(i)(t) 
%
%       where I(i,1) is a m x 1 vector with all zeros, except in position
%       (i,1), where I(i,1) = 1 (see also fOne). The Kronecker delta
%       delta(t,0) is one if t=0 and zero otherwise.
%
%	Output parameters:
%		out : p x N x nactive Jacobian matrix with the partial derivatives
%		      of y(t) w.r.t. the active elements in u0, where p is the
%		      number of outputs, N is the number of samples, and nactive is
%		      the number of active elements in u0.
%
%	Input parameters:
%       B_EdwxIdu0 : n x m matrix B + E partial(zeta(0))/partial(u0)
%		Edwdx : n x n x N matrix with N samples of A + E partial(zeta(t))/partial(x(t))
%		C : p x n matrix C
%		Fdwdx : p x n x N matrix with N samples of F partial(eta(t))/partial(x(t))
%		active : nactive x 1 vector with the indices in u0 indicating the
%		         active elements of u0 (active elements = elements
%		         on which optimization will be done).
%
%	Versions:
%       1.0 : April 20, 2016
%
%	Copyright (c) Vrije Universiteit Brussel – dept. ELEC
%   All rights reserved.
%   Software can be used freely for non-commercial applications only.
%   Disclaimer: This software is provided “as is” without any warranty.
%
%   See also fOne

%--------------------------------------------------------------------------
% Version 1.0
%--------------------------------------------------------------------------
% {

[p,~,N] = size(Fdwdx); % Number of outputs and number of samples
nactive = length(active); % Number of active parameters in u0

% Add C to F partial(eta(t))/partial(x(t)) for all samples at once (faster than doing this inside the loop)
if nactive > 0
    Fdwdx = Fdwdx + repmat(C,[1,1,N]);
end

out = zeros(p,N,nactive); % Preallocate
for k = 1:nactive % Loop over all active elements in u0
    Jprev = B_EdwxIdu0(:,active(k)); % Initialize partial derivative of x(1) w.r.t. u0(active(k))
    for t = 2:N % Loop over all samples
        % Calculate state update alternative state-space model at time t 
        J = Edwdx(:,:,t-1)*Jprev; % Terms in alternative states at time t-1
        % Calculate output alternative state-space model at time t
        out(:,t,k) = Fdwdx(:,:,t)*J;
        % Update previous state alternative state-space model
        Jprev = J;
    end
end

%}