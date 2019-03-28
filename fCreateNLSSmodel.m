function model = fCreateNLSSmodel(A,B,C,D,nx,ny,T1,T2,sat)
%FCREATENLSSMODEL Create polynomial nonlinear state-space model from initial linear state-space model.
%
%	Usage:
%		model = fCreateNLSSmodel(A,B,C,D,nx,ny,T1,T2,sat)
%
%	Description:
%		Create a polynomial nonlinear state-space model from a linear
%		initialization with state-space matrices A, B, C, and D. The state
%		equation is extended with a multivariate polynomial in the states
%		and the inputs. The nonlinear degree(s) of this polynomial is/are
%		specified in nx. Similarly, the output equation is extended with a
%		multivariate polynomial, where the degrees are specified in ny. The
%		transient handling is reflected in T1 (for periodic data) and T2
%		(for aperiodic data). A saturation nonlinearity instead of a
%		polynomial one is obsoleted; the optional parameter sat should be
%		zero if specified.
%
%	Output parameters:
%		model : structure containing the parameters and relevant data of
%		        the polynomial nonlinear state-space model. This structure
%		        has the following fields:
%               A : n x n state matrix
%               B : n x m input matrix
%               C : p x n output matrix
%               D : p x m feed-through matrix
%               lin.A : n x n state matrix of the linear initialization
%               lin.B : n x m input matrix of the linear initialization
%               lin.C : p x n output matrix of the linear initialization
%               lin.D : p x m feed-through matrix of the linear initialization
%               nx : vector with nonlinear degrees in state update
%               ny : vector with nonlinear degrees in output equation
%               n : number of states
%               m : number of inputs
%               p : number of outputs
%               xpowers : n_nx x (n+m) matrix containing the exponents of
%                         each of the n_nx monomials in the state update
%                         (see also fCombinations)
%               n_nx : number of monomials in the state update
%               E : n x n_nx matrix with polynomial coefficients in the
%                   state update
%               xactive : linear indices of the active elements in the
%                         transpose of the E matrix (active elements = 
%                         elements on which optimization will be done). By
%                         default, all elements in the E matrix are set as
%                         active. See fSelectActive to change this
%                         property.
%               ypowers : n_ny x (n+m) matrix containing the exponents of
%                         each of the n_ny monomials in the output equation
%                         (see also fCombinations)
%               n_ny : number of monomials in the output equation
%               F : p x n_ny matrix with polynomial coefficients in the
%                   output equation
%               yactive : linear indices of the active elements in the
%                         transpose of the F matrix (active elements = 
%                         elements on which optimization will be done). By
%                         default, all elements in the F matrix are set as
%                         active. See fSelectActive to change this
%                         property.
%               T1 : vector that indicates how the transient is handled for
%                    periodic signals (see also the Input parameters) 
%               T2 : scalar indicating how many samples from the start are
%                    removed or vector indicating which samples are removed
%                    (see also fComputeIndicesTransientRemovalArb)
%               sat : obsolete, zero flag
%               satCoeff : obsolete, n x 1 vector of ones
%
%	Input parameters:
%		A : n x n state matrix
%       B : n x m input matrix
%       C : p x n output matrix
%       D : p x m feed-through matrix
%       nx : vector with nonlinear degrees in state equation
%       ny : vector with nonlinear degrees in output equation
%       T1 : vector that indicates how the transient is handled for
%            periodic signals. The first element T1(1) is the number of
%            transient samples that should be prepended to each input
%            realization. The other elements T1(2:end) indicate the
%            starting sample of each realization in the signal. If T1 has
%            only one element, T1(2) is put to one (see also
%            fComputeIndicesTransient).
%       T2 : scalar indicating how many samples from the start are removed
%		     or vector indicating which samples are removed (see also
%		     fComputeIndicesTransientRemovalArb)
%       sat : obsolete, sat should be put equal to zero (optional, default
%             is zero)
%
%	Example:
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
%		1.0 : September 3, 2015
%       1.1 : April 20, 2016
%           Help updated
%       1.2 : August 12, 2016
%           Help updated
%
%	Copyright (c) Vrije Universiteit Brussel – dept. ELEC
%   All rights reserved.
%   Software can be used freely for non-commercial applications only.
%   Disclaimer: This software is provided “as is” without any warranty.
%
%	See also fComputeIndicesTransient, fComputeIndicesTransientRemovalArb, fFilterNLSS, fSScheckDims, fCombinations, fSelectActive

%--------------------------------------------------------------------------
% Version 1.0 = Version 1.1 = Version 1.2
%--------------------------------------------------------------------------
% {

% Linear state-space matrices
model.A = A;
model.B = B;
model.C = C;
model.D = D;

% State-space matrices linear intialization
model.lin.A = A;
model.lin.B = B;
model.lin.C = C;
model.lin.D = D;

% Nonlinear degrees
model.nx = nx; % in state equation
model.ny = ny; % in output equation

% State-space dimensions + checking their consistency
[model.n,model.m,model.p] = fSScheckDims(A,B,C,D);

% Nonlinear terms in state equation
model.xpowers = fCombinations(model.n+model.m,model.nx); % All possible terms
model.n_nx = size(model.xpowers,1); % Number of terms
model.E = zeros(model.n,model.n_nx); % Polynomial coefficients
model.xactive = (1:numel(model.E))'; % Active terms

% Nonlinear terms in state equation
model.ypowers = fCombinations(model.n+model.m,ny); % All possible terms
model.n_ny = size(model.ypowers,1); % Number of terms
model.F = zeros(model.p,model.n_ny); % Polynomial coefficients
model.yactive = (1:numel(model.F))'; % Active terms

% Transient handling
model.T1 = T1; % for periodic data
model.T2 = T2; % for aperiodic data

% Saturating nonlinearity or not (obsolete, sat should be zero)
if nargin > 8 && sat ~= 0
    warning('Saturating nonlinearity is obsolete, sat should be zero')
end
sat = 0;
model.sat = sat;
model.satCoeff = ones(model.n,1);

%}