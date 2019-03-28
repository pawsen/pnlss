function out = fJx0(A_Edwdx_0,Edwdx,C,Fdwdx,active)
%FJX0 Compute Jacobian w.r.t. x0 by filtering an alternative state-space model.
%
%	Usage:
%		out = fJx0(A_Edwdx_0,Edwdx,C,Fdwdx,active)
%
%	Description:
%		Computing the Jacobian of the output y(t) of a nonlinear
%		state-space model
%
%           x(t+1) = A x(t) + B u(t) + E zeta(x(t),u(t))
%           y(t)   = C x(t) + D u(t) + F eta(x(t),u(t))
%
%       w.r.t. the elements in the initial state x0 can be performed by
%       filtering an alternative nonlinear state-space model.
%       Let Jx0(i)(t) be the partial derivative of y(t) w.r.t. x0(i) and
%       let xx0(i)(t) be the partial derivative of x(t) w.r.t. x0(i).
%       Then the Jacobian Jx0 can be computed by filtering an alternative
%       state-space model
%
%           xx0(i)(t+1) = A xx0(i)(t) + E partial(zeta(t))/partial(x(t)) xx0(i)(t)
%           Jx0(i)(t)   = C xx0(i)(t) + F partial(eta(t))/partial(x(t))  xx0(i)(t)
%
%           with initial state xx0(i)(0) = I(i,1)
%
%       where I(i,1) is a n x 1 vector with all zeros, except in position
%       (i,1), where I(i,1) = 1 (see also fOne).
%
%	Output parameters:
%		out : p x N x nactive Jacobian matrix with the partial derivatives
%		      of y(t) w.r.t. the active elements in x0, where p is the
%		      number of outputs, N is the number of samples, and nactive is
%		      the number of active elements in x0.
%
%	Input parameters:
%       A_Edwdx_0 : n x n matrix A + E partial(zeta(0))/partial(x0)
%		Edwdx : n x n x N matrix with N samples of A + E partial(zeta(t))/partial(x(t))
%		C : p x n matrix C
%		Fdwdx : p x n x N matrix with N samples of F partial(eta(t))/partial(x(t))
%		active : nactive x 1 vector with the indices in x0 indicating the
%		         active elements of x0 (active elements = elements
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
nactive = length(active); % Number of active parameters in x0

% Add C to F partial(eta(t))/partial(x(t)) for all samples at once (faster than doing this inside the loop)
if nactive > 0
    Fdwdx = Fdwdx + repmat(C,[1,1,N]);
end

out = zeros(p,N,nactive); % Preallocate
for k = 1:nactive % Loop over all active elements in x0
    Jprev = A_Edwdx_0(:,active(k)); % Initialize partial derivative of x(1) w.r.t. x0(active(k))
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