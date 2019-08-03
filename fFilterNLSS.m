function [y,states] = fFilterNLSS(model,u,varargin)
%FFILTERNLSS Calculate the output and the states of a nonlinear state-space model with transient handling.
%
%	Usage:
%       [y,states] = fFilterNLSS(model,u,x0,u0) % implemented since version 1.1
%       [y,states] = fFilterNLSS(model,u,x0)
%       [y,states] = fFilterNLSS(model,u)
%       y = fFilterNLSS(model,u)
%
%   Description:
%       y = fFilterNLSS(model,u) calculates the output y of a nonlinear
%       state-space model by applying an input u and starting from a zero
%       initial state x0 and a zero initial input u0.
%       y = fFilterNLSS(model,u,x0) starts from an initial state x0.
%       y = fFilterNLSS(model,u,x0,u0) starts from an initial input u0
%       (implemented since version 1.1).
%       [y,states] = fFilterNLSS(model,u,x0, u0) also passes the states at 
%       each time step.
%
%   Output parameters:
%       y : N x p matrix with N samples of the p outputs
%       states : N x n matrix with N samples of the n states
%
%   Input parameters:
%       model : PNLSS model (see also fCreateNLSSmodel)
%       u : m x N matrix with N samples of the m inputs
%           If u is an N x m matrix, then u is transposed.
%       x0 : n x 1 vector with the initial states (optional, default is
%            zero initial states)
%       u0 : m x 1 vector with the initial inputs (optional, default is
%            zero initial inputs, implemented since version 1.1)
%
%   Example:
%       n = 3; % Number of states
%       m = 1; % Number of inputs
%       p = 1; % Number of outputs
%       sys = drss(n,p,m); % Random linear state-space model
%       nx = [2 3]; % Quadratic and cubic terms in state equation
%       ny = [2 3]; % Quadratic and cubic terms in output equation
%       T1 = 0; % No transient handling
%       T2 = []; % No transient handling
%       sat = 0; % Obsolete parameter sat = 0
%       model = fCreateNLSSmodel(sys.a,sys.b,sys.c,sys.d,nx,ny,T1,T2,sat); % Linear state-space model
%       N = 1e3; % Number of samples
%       u = randn(N,1); % Input signal
%       y = fFilterNLSS(model,u); % Modeled output signal
%       t = 0:N-1; % Time vector
%       y_lsim = lsim(sys,u,t); % Alternative way to calculate output of linear state-space model
%       figure
%           plot(t,y_lsim,'b')
%           hold on
%           plot(t,y,'r')
%           xlabel('Time')
%           ylabel('Output')
%           legend('lsim','PNLSS')
%
%	Versions:
%		1.0 : November 17, 2015
%       1.1 : December 15, 2015
%           Initial input can be specified
%           Debug: x0 = varargin{1} instead of varargin{3}
%       1.2 : April 19, 2016
%           Help updated
%
%	Copyright (c) Vrije Universiteit Brussel � dept. ELEC
%   All rights reserved.
%   Software can be used freely for non-commercial applications only.
%   Disclaimer: This software is provided �as is� without any warranty.
%
%	See also fCreateNLSSmodel, fComputeIndicesTransient, fFilterspeedNL, fComputeIndicesTransientRemoval

%--------------------------------------------------------------------------
% Version 1.1 = Version 1.2
%   Initial input can be specified
%   Debug: x0 = varargin{1} instead of varargin{3}
%--------------------------------------------------------------------------
% {

% Determine the initial state
if nargin >= 3 % If initial state specified
    x0 = varargin{1}; % Initial state
    if size(x0,1) ~= model.n || size(x0,2) ~= 1 % Check dimensions
        x0 = x0.';
    end
    if size(x0,1) ~= model.n || size(x0,2) ~= 1 % Check dimensions
        error('Initial state should have dimension n x 1')
    end
else % If initial state not specified
    x0 = zeros(model.n,1); % Use zero initial conditions
end

if nargin >= 4 % If initial input specified
    u0 = varargin{2}; % Initial input
    if size(u0,1) ~= model.m || size(u0,2) ~= 1 % Check dimensions
        error('Initial input should have dimension m x 1')
    end
else % If initial input not specified
    u0 = zeros(model.m,1); % Use zero initial input
end

% Check input dimensions
if size(u,1) ~= size(model.B,2) % If u has not as many rows as inputs
    if size(u,2) == size(model.B,2) % If u has as many columns as inputs
        % warning('Input matrix will be transposed to comply with the number of inputs in the model')
        u = u.'; % Transpose u
    else % If both the number of rows and columns of u are different from the number of inputs
        error('Number of inputs is inconsistent')
    end
end

% Collect model parameters and properties
[A,B,C,D,E,F]     = deal(model.A,model.B,model.C,model.D,model.E,model.F); % State-space matrices
[xpowers,ypowers] = deal(model.xpowers,model.ypowers); % All possible terms in state and output equation
[nx,ny]           = deal(model.nx,model.ny); % Nonlinear degrees in state and output equation
if nx == 0 % If no nonlinear terms in state equation
    nx = 1; % nx = 1 is equivalent to nx = 0
end
if ny == 0 % If no nonlinear terms in output equation
    ny = 1; % ny = 1 is equivalent to ny = 0
end

% Handle saturating nonlinearity (obsolete)
try
    sat = model.sat;
%     satCoeff = model.satCoeff;
catch
    sat = 0;
end

% Prepend transient samples to the input
N = size(u,2); % Number of samples
u = u(:,fComputeIndicesTransient(model.T1,N)); % Input with prepended transient samples

% Perform the actual filtering
if ~sat % If atan saturation is turned off
    [y,states] = fFilterspeedNL(A,B,C,D, ...
        E,F,xpowers,ypowers, ...
        max(nx),max(ny),u,x0,u0);
else % If saturation is on (obsolete)
%     [y,states] = filterspeedNLsatv2(A,B,C,D, ...
%         E,F, xpowers, ypowers, ...
%         max(nx),max(ny),satCoeff,u,x0);
    error('Saturation nonlinearity is obsolete')
end

% Transient handling
if model.T1 ~= 0 % If there were prepended transient samples
    states = states(:,fComputeIndicesTransientRemoval(model.T1,N,1)).'; % Remove transient samples from the states
    y = y(:,fComputeIndicesTransientRemoval(model.T1,N,1)).'; % Remove transient samples from the output
else % If no prepended transient samples
    states = states.'; % States (N x n)
    y = y.'; % Output (N x p)
end

%}

%--------------------------------------------------------------------------
% Version 1.0
%--------------------------------------------------------------------------
%{

% Determine the initial state
if nargin == 3 % If initial state specified
    x0 = varargin{3}; % Initial state
    if size(x0,1) ~= model.n || size(x0,2) ~= 1 % Check dimensions
        error('Initial state should have dimension n x 1')
    end
else % If initial state not specified
    x0 = zeros(model.n,1); % Use zero initial conditions
end

% Check input dimensions
if size(u,1) ~= size(model.B,2) % If u has not as many rows as inputs
    if size(u,2) == size(model.B,2) % If u has as many columns as inputs
        % warning('Input matrix will be transposed to comply with the number of inputs in the model')
        u = u.'; % Transpose u
    else % If both the number of rows and columns of u is different from the number of inputs
        error('Number of inputs is inconsistent')
    end
end

% Collect model parameters and properties
[A,B,C,D,E,F]     = deal(model.A,model.B,model.C,model.D,model.E,model.F); % State-space matrices
[xpowers,ypowers] = deal(model.xpowers,model.ypowers); % All possible terms in state and output equation
[nx,ny]           = deal(model.nx,model.ny); % Nonlinear degrees in state and output equation
if nx == 0 % If no nonlinear terms in state equation
    nx = 1; % nx = 1 is equivalent to nx = 0
end
if ny == 0 % If no nonlinear terms in output equation
    ny = 1; % ny = 1 is equivalent to ny = 0
end

% Handle saturating nonlinearity (obsolete)
try
    sat = model.sat;
%     satCoeff = model.satCoeff;
catch
    sat = 0;
end

% Prepend transient samples to the input
N = size(u,2); % Number of samples
u = u(:,fComputeIndicesTransient(model.T1,N)); % Input with prepended transient samples

% Perform the actual filtering
if ~sat % If atan saturation is turned off
    [y,states] = fFilterspeedNL(A,B,C,D, ...
        E,F,xpowers,ypowers, ...
        max(nx),max(ny),u,x0);
else % If saturation is on (obsolete)
%     [y,states] = filterspeedNLsatv2(A,B,C,D, ...
%         E,F, xpowers, ypowers, ...
%         max(nx),max(ny),satCoeff,u,x0);
    error('Saturation nonlinearity is obsolete')
end

% Transient handling
if model.T1 ~= 0 % If there were prepended transient samples
    states = states(:,fComputeIndicesTransientRemoval(model.T1,N,1)).'; % Remove transient samples from the states
    y = y(:,fComputeIndicesTransientRemoval(model.T1,N,1)).'; % Remove transient samples from the output
else % If no prepended transient samples
    states = states.'; % States (N x n)
    y = y.'; % Output (N x p)
end

%}
