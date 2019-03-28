%SJACOBIANANALYTICAL Script to compute the Jacobians in a nonlinear state-space model.
%
%	Description:
%       sJacobianAnalytical is a script that computes the steady-state
%       Jacobians of a nonlinear state-space model
%
%           x(t+1) = A x(t) + B u(t) + E zeta(x(t),u(t))
%           y(t)   = C x(t) + D u(t) + F eta(x(t),u(t))
%
%       i.e. the partial derivatives of the modeled output w.r.t. the
%       active elements in the A, B, E, F, D, and C matrices.
%       This script is called in fLMnlssWeighted.
%
%	Versions:
%		1.0 : November 16, 2015
%       1.1 : April 20, 2016
%           Debug: commented out the check for active elements in the E and
%           F matrix (the empty zeta matrix would otherwise produce an
%           error if there are no active elements in E).
%           Help updated.
%
%	Copyright (c) Vrije Universiteit Brussel – dept. ELEC
%   All rights reserved.
%   Software can be used freely for non-commercial applications only.
%   Disclaimer: This software is provided “as is” without any warranty.
%
%   See also fLMnlssWeighted, fComputeIndicesTransient, fEdwdx, fTermNL, fComputeJF, fComputeIndicesTransientRemoval, fJNL

%--------------------------------------------------------------------------
% Version 1.1
%   Debug: commented out the check for active elements in the E and F
%   matrix (the empty zeta matrix would otherwise produce an error if there
%   are no active elements in E).
%   Help updated.
%--------------------------------------------------------------------------
% {

% Prepend transient samples to the states
statesT = states(fComputeIndicesTransient(T1,N),:);

% Compute A + E partial(zeta(t))/partial(x(t)) and zeta
% if ~isempty(xactive) % If active elements in E matrix
    contrib = [statesT uT].'; % Collect states and outputs with prepended transient samples
    degx = max(max(xpowers)); % Maximum degree in individual state or input in zeta
    A_EdwxIdx = fEdwdx(contrib,xd_powers,squeeze(xd_coeff),E,degx,n) ...
                + repmat(A,[1 1 NT]); % E partial(zeta(t))/partial(x(t)) + A (n x n x NT)
    zeta = fTermNL([statesT uT].',xpowers,degx).'; % zeta (NT x n_nx)
% else % If no active elements in E matrix
%     zeta = []; % No zeta
%     A_EdwxIdx = repmat(A,[1 1 NT]); % Zero E matrix => only A
% end

% Compute F partial(eta(t))/partial(x(t)), eta, and Jacobian w.r.t. F
% if ~isempty(yactive) % If active elements in F matrix
    contrib = [statesT uT].'; % Collect states and outputs with prepended transient samples
    degy = max(max(ypowers)); % Maximum degree in individual state or input in eta
    FdwyIdx = fEdwdx(contrib,yd_powers,squeeze(yd_coeff),F,degy,n); % F partial(eta(t))/partial(x(t)) (p x n x NT)
    eta = fTermNL(contrib,ypowers,degy); % eta (n_ny x NT)
    JF = fComputeJF(p,yactive,n_ny,eta); % Jacobian w.r.t. F with prepended transient samples (p*NT x nactiveF)
    JF = JF(fComputeIndicesTransientRemoval(T1,N,p),:); % Remove transient samples (p*N x nactiveF)
    clear eta; % Free some memory
% else % If no active elements in F matrix
%     FdwyIdx = zeros(p,n,NT); % Zero F matrix
%     JF = []; % No Jacobian w.r.t. F
% end

% Compute Jacobian w.r.t. E
JE = fJNL(zeta,A_EdwxIdx,C,FdwyIdx,xactive(:)); % p x NT x nactive
JE = reshape(permute(JE,[2 1 3]),p*NT,length(xactive)); % p*NT x nactive
JE = JE(fComputeIndicesTransientRemoval(T1,N,p),:); % Remove transient samples (p*N x nactive)
clear zeta; % zeta no longer needed

% Compute Jacobian w.r.t. A
JA = fJNL(statesT,A_EdwxIdx,C,FdwyIdx,(1:n^2).'); % p x NT x n^2
JA = reshape(permute(JA,[2 1 3]),p*NT,n^2); % p*NT x n^2
JA = JA(fComputeIndicesTransientRemoval(T1,N,p),:); % Remove transient samples (p*N x n)

% Compute Jacobian w.r.t. B
JB = fJNL(uT,A_EdwxIdx,C,FdwyIdx,(1:n*m).'); % p x NT x n*m
JB = reshape(permute(JB,[2 1 3]),p*NT,n*m); % p*NT x n*m
JB = JB(fComputeIndicesTransientRemoval(T1,N,p),:); % Remove transient samples (p*N x n*m)

% Compute transient-free Jacobian w.r.t. active elements in A, B, E, F, D, and C
JD = kron(eye(p),u); % Jacobian w.r.t. D (p*N x p*m)
JC = kron(eye(p),states); % Jacobian w.r.t. C (p*N x p*n)
J = [JA JB JE JF JD JC]; % Put together all Jacobians
J = J(without_T2,:); % Remove transient points from nonperiodic transient handling
clear A_EdwxIdx FdwyIdx JA JB JC JD JE JF; % Free some memory

%}

%--------------------------------------------------------------------------
% Version 1.0
%--------------------------------------------------------------------------
%{

% Prepend transient samples to the states
statesT = states(fComputeIndicesTransient(T1,N),:);

% Compute A + E partial(zeta(t))/partial(x(t)) and zeta
if ~isempty(xactive) % If active elements in E matrix
    contrib = [statesT uT].'; % Collect states and outputs with prepended transient samples
    degx = max(max(xpowers)); % Maximum degree in individual state or input in zeta
    A_EdwxIdx = fEdwdx(contrib,xd_powers,squeeze(xd_coeff),E,degx,n) ...
                + repmat(A,[1 1 NT]); % E partial(zeta(t))/partial(x(t)) + A (n x n x NT)
    zeta = fTermNL([statesT uT].',xpowers,degx).'; % zeta (NT x n_nx)
else % If no active elements in E matrix
    zeta = []; % No zeta
    A_EdwxIdx = repmat(A,[1 1 NT]); % Zero E matrix => only A
end

% Compute F partial(eta(t))/partial(x(t)), eta, and Jacobian w.r.t. F
if ~isempty(yactive) % If active elements in F matrix
    contrib = [statesT uT].'; % Collect states and outputs with prepended transient samples
    degy = max(max(ypowers)); % Maximum degree in individual state or input in eta
    FdwyIdx = fEdwdx(contrib,yd_powers,squeeze(yd_coeff),F,degy,n); % F partial(eta(t))/partial(x(t)) (p x n x NT)
    eta = fTermNL(contrib,ypowers,degy); % eta (n_ny x NT)
    JF = fComputeJF(p,yactive,n_ny,eta); % Jacobian w.r.t. F with prepended transient samples (p*NT x nactiveF)
    JF = JF(fComputeIndicesTransientRemoval(T1,N,p),:); % Remove transient samples (p*N x nactiveF)
    clear eta; % Free some memory
else % If no active elements in F matrix
    FdwyIdx = zeros(p,n,NT); % Zero F matrix
    JF = []; % No Jacobian w.r.t. F
end

% Compute Jacobian w.r.t. E
JE = fJNL(zeta,A_EdwxIdx,C,FdwyIdx,xactive(:)); % p x NT x nactive
JE = reshape(permute(JE,[2 1 3]),p*NT,length(xactive)); % p*NT x nactive
JE = JE(fComputeIndicesTransientRemoval(T1,N,p),:); % Remove transient samples (p*N x nactive)
clear zeta; % zeta no longer needed

% Compute Jacobian w.r.t. A
JA = fJNL(statesT,A_EdwxIdx,C,FdwyIdx,(1:n^2).'); % p x NT x n^2
JA = reshape(permute(JA,[2 1 3]),p*NT,n^2); % p*NT x n^2
JA = JA(fComputeIndicesTransientRemoval(T1,N,p),:); % Remove transient samples (p*N x n)

% Compute Jacobian w.r.t. B
JB = fJNL(uT,A_EdwxIdx,C,FdwyIdx,(1:n*m).'); % p x NT x n*m
JB = reshape(permute(JB,[2 1 3]),p*NT,n*m); % p*NT x n*m
JB = JB(fComputeIndicesTransientRemoval(T1,N,p),:); % Remove transient samples (p*N x n*m)

% Compute transient-free Jacobian w.r.t. active elements in A, B, E, F, D, and C
JD = kron(eye(p),u); % Jacobian w.r.t. D (p*N x p*m)
JC = kron(eye(p),states); % Jacobian w.r.t. C (p*N x p*n)
J = [JA JB JE JF JD JC]; % Put together all Jacobians
J = J(without_T2,:); % Remove transient points from nonperiodic transient handling
clear A_EdwxIdx FdwyIdx JA JB JC JD JE JF; % Free some memory

%}