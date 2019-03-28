function [y,states] = fFilterspeedNL(A,B,C,D,E,F,xpowers,ypowers,max_nx,max_ny,u,x0,u0)
%FFILTERSPEEDNL Calculate the output and the states of a nonlinear state-space model without transient handling.
%
%	Usage:
%		[y,states] = fFilterspeedNL(A,B,C,D,E,F,xpowers,ypowers,max_nx,max_ny,u,x0,u0) % implemented since version 1.1
%		[y,states] = fFilterspeedNL(A,B,C,D,E,F,xpowers,ypowers,max_nx,max_ny,u,x0)
%
%	Description:
%		Calculate the output and the states of a nonlinear state-space
%		model
%           x(t+1) = A x(t) + B u(t) + E zeta(x(t),u(t))
%           y(t)   = C x(t) + D u(t) + F eta(x(t),u(t))
%       where zeta and eta are polynomials whose exponents are given in
%       xpowers and ypowers, respectively. The maximum degree in one
%       variable (a state or an input) in zeta or eta is given in max_nx
%       and max_ny, respectively. The initial state is given in x0.
%       This is a low-level function that is used in fFilterNLSS.
%
%	Output parameters:
%       y : p x N matrix with N samples of the p outputs
%       states : n x N matrix with N samples of the n states
%
%	Input parameters:
%		A : n x n state matrix
%       B : n x m input matrix
%       C : p x n output matrix
%       D : p x m feed-through matrix
%		E : n x nzeta matrix with polynomial coefficients in the state
%		    equation
%		F : p x neta matrix with polynomial coefficients in the output
%		    equation 
%		xpowers : nzeta x (n+m) matrix with the exponents in the states and
%		          the inputs for each of the nzeta terms
%		ypowers : neta x (n+m) matrix with the exponents in the states and
%		          the inputs for each of the neta terms 
%		max_nx : maximum degree in one variable in zeta
%		max_ny : maximum degree in one variable in eta
%		u : m x N matrix with N samples of the m inputs
%		x0 : n x 1 vector with the initial state
%       u0 : m x 1 vector with the initial input (optional, default is zero
%            initial inputs, implemented since version 1.1)
%
%   Example:
%       % A PNLSS model with only input nonlinearities has an equivalent
%       % linear state-space model structure with an extended input vector.
%       % Output of a PNLSS model with only input nonlinearities
%       sys = drss(2,1,1); % Second-order SISO discrete-time state-space model
%       set(sys,'Ts',1); % Set unit sampling time
%       [A,B,C,D] = ssdata(sys); % Linear state-space matrices
%       xpowers = [0 0 2;
%                  0 0 3]; % Quadratic and cubic input term in state update
%       ypowers = [0 0 2;
%                  0 0 3]; % Quadratic and cubic input term in output equation
%       E = [1 2;
%            3 4]; % Monomial coefficients in state update
%       F = [1 2]; % Monomial coefficients in output equation
%       max_nx = 3; % Maximum degree of a monomial in state update
%       max_ny = 3; % Maximum degree of a monomial in output equation
%       u = randn(1,1000); % Input signal
%       x0 = [1; 2]; % Start from non-zero initial state
%       y = fFilterspeedNL(A,B,C,D,E,F,xpowers,ypowers,max_nx,max_ny,u,x0); % Compute output
%       % Output corresponding linear state-space model with extended input vector
%       u_ext = [u(:) u(:).^2 u(:).^3]; % Extended input vector
%       B_ext = [B E]; % Extended input matrix
%       D_ext = [D F]; % Extended feed-through matrix
%       t = 0:999; % Time vector
%       y_lsim = lsim(ss(A,B_ext,C,D_ext,1),u_ext,t,x0); % Output extended linear state-space model
%       % Compare the PNLSS and lsim approach to calculate the output
%       figure
%           plot(t,y_lsim,'b')
%           hold on
%           plot(t,y,'r')
%           xlabel('Time')
%           ylabel('Output')
%           legend('lsim','PNLSS')
%
%	Versions:
%		1.0 : November 27, 2015
%       1.1 : December 15, 2015
%           Initial input can be specified
%       1.2 : April 19, 2016
%           Help updated
%
%	Copyright (c) Vrije Universiteit Brussel – dept. ELEC
%   All rights reserved.
%   Software can be used freely for non-commercial applications only.
%   Disclaimer: This software is provided “as is” without any warranty.
%
%	See also fFilterNLSS

%--------------------------------------------------------------------------
% Version 1.1 = Version 1.2
%   Initial input can be specified
%--------------------------------------------------------------------------
% {

% Set default value for initial input if not specified
if (nargin < 13) || isempty(u0)
    u0 = zeros(size(u,1),1); % Zero initial input
end

x_t = x0; % Initial state
u_t = u0; % Initial input
y = NaN(size(C,1),size(u,2)); % Preallocate
states = zeros(size(A,1),size(u,2)); % Preallocate
repmat_x = ones(size(xpowers,1),1); % Faster than repmat in the loop
repmat_y = ones(size(ypowers,1),1); % Faster than repmat in the loop
for t = 1:size(u,2) % Loop over all samples
    % State equation x(t) = A*x(t-1) + B*u(t-1) + E*zeta(x(t-1),u(t-1))
    zeta_t = prod((repmat_x*[x_t; u_t].').^xpowers,2);
    x_t = [A B E]*[x_t; u_t; zeta_t];
    
    % Output equation y(t) = C*x(t) + D*u(t) + F*eta(x(t),u(t))
    u_t = u(:,t);
    eta_t = prod((repmat_y*[x_t; u_t].').^ypowers,2);
    y(:,t) = [C D F]*[x_t; u_t; eta_t];
    
    % Collect states
    states(:,t) = x_t;
end

%}

%--------------------------------------------------------------------------
% Version 1.0
%--------------------------------------------------------------------------
%{

x_t = x0; % Initial state
u_t = zeros(size(u,1),1); % Initial input
y = NaN(size(C,1),size(u,2)); % Preallocate
states = zeros(size(A,1),size(u,2)); % Preallocate
repmat_x = ones(size(xpowers,1),1); % Faster than repmat in the loop
repmat_y = ones(size(ypowers,1),1); % Faster than repmat in the loop
for t = 1:size(u,2) % Loop over all samples
    % State equation x(t) = A*x(t-1) + B*u(t-1) + E*zeta(x(t-1),u(t-1))
    zeta_t = prod((repmat_x*[x_t; u_t].').^xpowers,2);
    x_t = [A B E]*[x_t; u_t; zeta_t];
    
    % Output equation y(t) = C*x(t) + D*u(t) + F*eta(x(t),u(t))
    u_t = u(:,t);
    eta_t = prod((repmat_y*[x_t; u_t].').^ypowers,2);
    y(:,t) = [C D F]*[x_t; u_t; eta_t];
    
    % Collect states
    states(:,t) = x_t;
end

%}