function [model,y_mod,models,Cost] = fLMnlssWeighted_x0u0(u,y,model,MaxCount,W,lambda,LambdaJump)
%FLMNLSSWEIGHTED_X0U0 Optimize PNLSS model and initial conditions using weighted Levenberg-Marquardt algorithm.
%
%	Usage:
%		[model,y_mod,models,Cost] = fLMnlssWeighted_x0u0(u,y,model,MaxCount,W,lambda,LambdaJump)
%       [model,y_mod,models,Cost] = fLMnlssWeighted_x0u0(u,y,model,MaxCount,W)
%
%	Description:
%		fLMnlssWeighted_x0u0 performs a Levenberg-Marquardt optimization on
%		the parameters of a PNLSS model (i.e. the elements of the matrices
%		A, B, C, D, E, and F) and on the initial conditions (i.e. the
%		initial state x0 = x(0) and the initial input u0 = u(0)). The
%		difference between the modeled and the measured output is minimized
%		in a weighted least squares sense, either in the time or the
%		frequency domain. A simple stabilization method can be applied by
%		simulating a validation data set during estimation and checking
%		whether or not the modeled output stays within prespecified bounds.
%		If not, the Levenberg-Marquardt iteration acts as if the cost
%		function increased.
%
%	Output parameters:
%		model : optimized model (= best on estimation data)
%       y_mod : output of the optimized model
%       models : collection of models (initial model + model after a
%                successful iteration)
%       Cost : collection of the unweighted rms error at each iteration
%              (NaN if iteration was not successful (i.e. when the weighted
%              rms error increased))
%
%	Input parameters:
%		u : N x m input signal
%       y : N x p output signal
%       model : initial model (see fCreateNLSSmodel) with extra fields
%               x0 : initial value of the initial state (optional, default
%                    = zeros(n,1))
%               u0 : initial value of the initial input (optional, default
%                    = zeros(m,1))
%               x0active : vector with the indices of the active elements
%                         (i.e. elements that will be optimized) in x0
%                         (optional, default = [])
%               u0active : vector with the indices of the active elements
%                         (i.e. elements that will be optimized) in u0
%                         (optional, default = [])
%               Additionally, two extra optional fields can be added to the
%               model to perform a simple stabilization method (both fields
%               are needed to perform this method):
%               u_val : m x N_val matrix with N_val samples of the
%                       validation input (optional, no default value)
%               max_out : bound on the maximum absolute value of the
%                         simulated validation output (optional, no default
%                         value)
%               After each Levenberg-Marquardt iteration, the validation
%               data is simulated (without taking into account the
%               transient settings). If the simulated validation output
%               does not respect the max_out bound, then the iteration is
%               considered unsuccessful.
%       MaxCount : (maximum) number of iterations (there is not yet an
%                  early stopping criterion)
%       W : p x p x NFD weighting matrix if frequency-domain weighting
%           (e.g. square root of covariance matrix of output noise
%           spectrum), where NFD is the number of frequency bins in the
%           positive half of the spectrum (of one period and one phase
%           realization) of the input (e.g. NFD = floor(Npp/2), where Npp
%           is the number of samples in one period and one phase
%           realization for a multisine excitation).
%           N x p weighting sequences if time-domain weighting.
%           [] if no weighting.
%           (optional, default is no weighting)
%       lambda : initial Levenberg-Marquardt parameter
%                (optional, default = 0, which corresponds to a
%                Gauss-Newton algorithm). After a successful iteration,
%                lambda is halved. After an unsuccessful iteration, lambda
%                is multiplied with a factor sqrt(10), unless lambda was
%                zero, in which case lambda is put equal to the dominant
%                singular value of the Jacobian.
%       LambdaJump : each LambdaJump iterations, the Levenberg-Marquardt
%                    parameter is made smaller by a factor 10, so that the
%                    algorithm leans more towards a Gauss-Newton algorithm,
%                    which converges faster than a gradient-descent
%                    algorithm (optional, default = 1001)
%
%	Example:
%       % Model input/output data of a Hammerstein system with non-zero initial conditions 
% 		N = 200; % Number of samples
%       NTrans = 100; % Number of samples after zero initial conditions
%       u = randn(NTrans+N,1); % Input signal
%       f_NL = @(x) x + 0.2*x.^2 + 0.1*x.^3; % Nonlinear function
%       [b,a] = cheby1(2,5,2*0.3); % Filter coefficients
%       x = f_NL(u); % Intermediate signal
%       y = filter(b,a,x); % Output signal
%       u(1:NTrans) = []; % Remove first NTrans samples to obtain non-zero intial conditions 
%       x(1:NTrans) = []; % Remove first NTrans samples to obtain non-zero intial conditions
%       y(1:NTrans) = []; % Remove first NTrans samples to obtain non-zero intial conditions
%       scale = u\x; % Scale factor
%       sys = ss(tf(scale*b,a,[])); % Initial linear model = scale factor times underlying dynamics 
%       nx = [2 3]; % Quadratic and cubic terms in state equation
%       ny = [2 3]; % Quadratic and cubic terms in output equation
%       T1 = 0; % No periodic signal transient handling
%       T2 = []; % No transient samples to discard
%       model = fCreateNLSSmodel(sys.a,sys.b,sys.c,sys.d,nx,ny,T1,T2); % Initial linear model
%       model.xactive = fSelectActive('inputsonly',2,1,2,nx); % A Hammerstein system only has nonlinear terms in the input
%       model.yactive = fSelectActive('inputsonly',2,1,1,nx); % A Hammerstein system only has nonlinear terms in the input
%       MaxCount = 50; % Maximum number of iterations
%       W = []; % No weighting
%       [modelOpt,yOpt] = fLMnlssWeighted(u,y,model,MaxCount,W); % Optimized model and modeled output (without estimating initial conditions)
%       model_x0u0 = model; % Estimate initial conditions
%       model_x0u0.x0active = (1:model_x0u0.n).'; % Estimate initial conditions
%       model_x0u0.u0active = (1:model_x0u0.m).'; % Estimate initial conditions
%       [modelOpt_x0u0,yOpt_x0u0] = fLMnlssWeighted_x0u0(u,y,model_x0u0,MaxCount,W); % Optimized model and modeled output (initial conditions estimated);
%       t = 0:N-1;
%       figure
%           plot(t,y,'b')
%           hold on
%           plot(t,y-yOpt,'r')
%           plot(t,y-yOpt_x0u0,'g')
%           xlabel('Time')
%           ylabel('Output / output error')
%           legend('True','Error PNLSS','Error PNLSS (initial conditions estimated)')
%
%	Reference:
%		Paduart, J., Lauwers, L., Swevers, J., Smolders, K., Schoukens, J.,
%		and Pintelon, R. (2010). Identification of nonlinear systems using
%		Polynomial Nonlinear State Space models. Automatica, 46:647-656.
%
%	Versions:
%		1.0 : January 13, 2016
%       1.1 : April 20, 2016
%           Help updated
%       1.2 : September 1, 2016
%           Fixed a bug where the time-domain weighting produced a
%           matrix-valued cost function instead of a scalar one.
%
%	Copyright (c) Vrije Universiteit Brussel – dept. ELEC
%   All rights reserved.
%   Software can be used freely for non-commercial applications only.
%   Disclaimer: This software is provided “as is” without any warranty.
%
%	See also fCreateNLSSmodel, fFilterNLSS, fComputeIndicesTransientRemovalArb, fWeightJacobSubSpace, fReIm 

%--------------------------------------------------------------------------
% Version 1.2
%   Fixed a bug where the time-domain weighting produced a matrix-valued
%   cost function instead of a scalar one.
%--------------------------------------------------------------------------
% {

% Set default value for weighting if not specified
if nargin < 5
    W = []; % No weighting
end

% Set default value for Levenberg-Marquardt parameter if not specified
if (nargin < 6) || isempty(lambda)
    lambda = 0; % Later on calculated as the dominant singular value of the Jacobian
end

% Set default value for LambdaJump
% Each LambdaJump iterations, the Levenberg-Marquardt parameter is made
% smaller, so that the algorithm leans more towards a Gauss-Newton
% algorithm (converges faster) and less to a gradient descent method
% (converges in a larger range)
if (nargin < 7) || isempty(LambdaJump)
    LambdaJump = 1001;
end

% Extract parameters for later use from input arguments
N                 = length(u); % Number of samples
[A,C,E,F]         = deal(model.A,model.C,model.E,model.F); %#ok State-space matrices used in analytical calculation Jacobian
[n,m,p]           = deal(model.n,model.m,model.p); %#ok Number of states, inputs, and outputs used in analytical calculation Jacobian
[n_nx,n_ny]       = deal(model.n_nx,model.n_ny); % Number of terms in state and output equation
[xpowers,ypowers] = deal(model.xpowers,model.ypowers); % All possible terms in state and output equation
xactive           = model.xactive; %#ok Active terms in state equation used in analytical calculation Jacobian
yactive           = model.yactive; %#ok Active terms in output equation used in analytical calculation Jacobian
[T1,T2,sat]       = deal(model.T1, model.T2,model.sat); % Transient and saturation parameters

% Extract parameters initial conditions
B = model.B; %#ok Input matrix used in analytical calculation Jacobian
if isfield(model,'x0')
    x0 = model.x0; % Initial state vector
else
    x0 = zeros(n,1); % Zero initial state vector if not specified
    model.x0 = x0;
end
if isfield(model,'u0')
    u0 = model.u0; % Initial input
else
    u0 = zeros(m,1); % Zero initial input if not specified
    model.u0 = u0;
end
if isfield(model,'x0active')
    x0active = model.x0active; %#ok: Indices active elements in initial state vector used in analytical calculation Jacobian
else
    x0active = []; % By default no optimization on initial state vector
    model.x0active = x0active;
end
if isfield(model,'u0active')
    u0active = model.u0active; %#ok: Indices active elements in initial input vector used in analytical calculation Jacobian
else
    u0active = []; % By default no optimization on initial input vector
    model.u0active = u0active;
end

% Compute the (transient-free) modeled output and the corresponding states
[y_mod,states] = fFilterNLSS(model,u,x0,u0); %#ok: states used in analytical Jacobian calculation script

% Determine if weighting is in frequency or time domain: only implemented for periodic signals.
NFD = size(W,3); % Number of frequency bins where weighting is specified (e.g. NFD = floor(Npp/2), where Npp is the number of samples in one period and one phase realization for a multisine excitation)
if isempty(W)
    % No weighting
    FreqWeighting = false; % No frequency weighting
    W = ones(N,p); % No weighting
elseif NFD > 1
    % Frequency weighting
    FreqWeighting = true; % Frequency weighting
    R = round((N-T2)/NFD/2); % Number of realizations
    if mod(N-T2,R) ~= 0
        error('Transient handling and weighting matrix are incompatible')
    end
else
    % Time-domain weighting
    FreqWeighting = false; % No frequency weighting
end

% If T2 is a scalar, it denotes the number of transient points to discard.
% If it is a vector, it denotes the indices of the points to discard, e.g.
% when several data sequences were put together.
without_T2 = fComputeIndicesTransientRemovalArb(T2,N,p); % Samples that are not discarded

% Compute the (weighted) error signal
err_old = y_mod(without_T2) - y(without_T2); % Error signal (without transient for aperiodic signals)
if FreqWeighting
    err_old = reshape(err_old,[(N-T2)/R,R,p]);
    err_old = fft(err_old); % Go to the frequency domain
    err_old = err_old(1:NFD,:,:); % Select only the positive half of the spectrum
    err_old = permute(err_old,[3 1 2]); % p x NFD x R
    err_old = fWeightJacobSubSpace(err_old,W,p,NFD,R); % Add weighting
    K_old = fVec(err_old)'*fVec(err_old); % Calculate cost function
    err_old = permute(err_old,[2 3 1]); % NFD x R x p
    err_old = err_old(:); % NFD R p x 1
    err_old = fReIm(err_old); % Split in real and imaginary part
else
    err_old = err_old.*W(without_T2); % Time-domain weighting
    K_old   = fVec(err_old)'*fVec(err_old); % Compute cost function
end

% Initialization of Levenberg-Marquardt algorithm
Count = 1; % Iteration number
models = []; % Collection of all models after a successful step
Cost = NaN(MaxCount,1); % Sequence of cost functions after a successful step

% Compute the rms value of the weighted output to later on calculate the
% relative error after each successful iteration
if ~FreqWeighting
    rms_y = rms(W(without_T2).*y(without_T2));
end

% Compute the derivatives of the polynomials zeta and eta
[xd_powers,xd_coeff] = lfPolyDerivative(xpowers); %#ok Derivatives of zeta used in analytical calculation Jacobian
[yd_powers,yd_coeff] = lfPolyDerivative(ypowers); %#ok Derivatives of eta used in analytical calculation Jacobian

% Extract the transient part of the input
uT = u(fComputeIndicesTransient(T1,N),:); % Transient part of the input
NT = length(uT); % Length of the transient part

% Should the model be stable on the validation set?
if all(isfield(model,{'u_val','max_out'}))
    stabilisation = true; % Yes, model stable on validation set
else
    stabilisation = false; % No
end

% Prepare for Levenberg-Marquardt optimization
cwd = pwd; % The current working directory
cd(tempdir); % Switch to the system's temporary folder
warning('off','MATLAB:pack:InvalidInvocationLocation'); % Disable warning
pack; % Consolidate workspace memory
disp('Starting L.M. Optimization...')
% inittime = clock; % Save current time to estimate the end time during the optimization 

% Levenberg-Marquardt optimization
while Count < MaxCount
    % Compute the Jacobian w.r.t. the elements in the matrices A, B, E, F, D, and C
    sJacobianAnalytical_x0u0;
    
    % Add frequency- or time-domain weighting
    if FreqWeighting
        numpar = size(J,2); % Number of parameters
        J = reshape(J,[(N-T2)/R,R,p,numpar]); % Npp x R x p x numpar
        J = fft(J); % Go to the frequency domain
        J = J(1:NFD,:,:,:); % Select only the positive half of the spectrum
        J = permute(J,[3 1 2 4]); % p x NFD x R x numpar
        J = reshape(J,[p NFD R*numpar]); % p x NFD x R numpar
        J = fWeightJacobSubSpace(J,W,p,NFD,R*numpar); % Add weighting
        J = reshape(J,[p NFD R numpar]); % p x NFD x R x numpar
        J = permute(J,[2 3 1 4]); % NFD x R x p x numpar
        J = reshape(J,[NFD*R*p,numpar]); % NFD R p x numpar
        J = fReIm(J); % Split in real and imaginary part
    else
        J = J.*repmat(W(without_T2),1,size(J,2)); % Add time-domain weighting
    end

    % (Potentially) improve numerical conditioning
    [J,scaling] = fNormalizeColumns(J); % Normalize columns to unit rms value

    K = K_old; % Initial value of the cost function (=> unsuccessful step)
    [U,S,V] = svd(J,0); % Singular value decomposition of the Jacobian
    clear J;
    pack; % Consolidate workspace memory

    if isnan(K)
        break % Stop nonlinear optimization (typically when initial model is unstable)
    end
    
    nowexit = false; % nowexit is true when a successfull step was found
    
    % Keep looking for a successful step in the parameter update
    while not(nowexit) && Count <= MaxCount % K >= K_old && Count <= MaxCount
        % Estimate the rank of the Jacobian
        s = diag(S); % Singular values of the Jacobian
        tol = NT*eps(max(s)); % Tolerance to determine the number of rank deficiencies
        r = min(sum(s>tol),length(s)-n^2); % Estimate the rank of the Jacobian (there are at least n^2 rank deficiencies)
        
        % Compute parameter update dtheta = -L*err_old, where L = (J.'*J + lambda*eye)^-1*J.' 
        s = diag(s(1:r)./(s(1:r).^2 + lambda^2*ones(r,1))); % Singular values of L
        dtheta = -V(:,1:r)*s*U(:,1:r).'*err_old; % Parameter update
        dtheta = dtheta./scaling'; % Denormalize the estimates (Jacobian was normalized)
        
        % Update the model and save as model_new
        model_new = lfVec2Par(dtheta,model);
        
        % Compute new model output and states
        [y_mod, states_new] = fFilterNLSS(model_new,u,model_new.x0,model_new.u0);
        
        % If the model should be stable on the validation set, ...
        if stabilisation
            % ... compute modeled output on validation set
            model_new.T1 = 0; % No transient handling on validation set
            yStab = fFilterNLSS(model_new,model.u_val,model_new.x0,model_new.u0); % Modeled output
            model_new.T1 = T1; % Transient handling on estimation set
        end
        % If model is not stable on validation set, while this was required, ... 
        if stabilisation && (max(abs(yStab(:))) > model.max_out || any(isnan(yStab(:)))>0)
            K = Inf; % ... set cost function to infinity (=> unsuccessful step), ...
            disp('Unstable for Validation Set');
        else
            % ... else, compute new model error
            err = y_mod(without_T2) - y(without_T2); % Discard transient
            rms_e = rms(err(:)); % Unweighted rms error
            if FreqWeighting
                err = reshape(err,[(N-T2)/R,R,p]); % Npp x R x p
                err = fft(err); % Go to the frequency domain
                err = err(1:NFD,:,:); % Select only the positive half of the spectrum
                err = permute(err,[3 1 2]); % p x Npp x R
                err = fWeightJacobSubSpace(err,W,p,NFD,R); % Add weighting
                K = fVec(err)'*fVec(err); % Weighted least-squares cost function
                err = permute(err,[2 3 1]); % Npp x R x p
                err = err(:); % Npp R p x 1
                err = fReIm(err); % Split in real and imaginary part
            else
                err = err.*W(without_T2); % Time domain weighting (Multiple Output: vec(err))
                K = fVec(err)'*fVec(err); % Cost function
            end
        end
        
        % Check successfullness parameter update
        if K >= K_old % Step was not succesful
            if lambda == 0 % lambda was initialed as zero (lambda = 0 should not be possible later on)
                lambda = S(1,1); % Initial Levenberg-Marquardt parameter = dominant singular value
            else
                lambda = lambda*sqrt(10); % Lean more towards gradient descent method (converges in larger range)
            end;
        elseif isnan(K) % Typically if model is unstable
            lambda = lambda*sqrt(10); % Lean more towards gradient descent method (converges in larger range)
        else % Step was succesful
            lambda = lambda/2; % Lean more towards Gauss-Newton algorithm (converges faster)
            nowexit = true; % Stop searching for successful step
            if FreqWeighting
                % Display weighted mean square cost function (divided by 2*NFD*R*p) and the condition number of the Jacobian 
                disp([num2str(Count) ' Cost Function: ' num2str(K/NFD/R/p/2) ' - Cond: ' num2str(S(1,1)/S(r,r))]);
            else
                % Display rms error normalized with weighted output, the condition number of the Jacobian, the Levenberg Marquardt parameter, and the unweighted rms error
                disp([num2str(Count) '  RelErr ' num2str(rms_e/rms_y) ' - Cond: ' num2str(S(1,1)/S(r,r)) ...
                    ' - Lambda: ' num2str(lambda) '  RmsErr ' num2str(rms_e)]);
            end
            Cost(Count) = rms_e; % Collect unweighted rms errors after each successful step
        end
        
        % Make Levenberg-Marquardt parameter smaller after LambdaJump iterations
        if (rem(Count,LambdaJump)==0 && Count~=0)
            lambda = lambda/10; % Lean more towards Gauss-Newton method (converges faster)
            disp('Forced jump lambda');
        end
        
        % Update the number of iterations
        Count = Count + 1;
    end
    clear U V;
    
    % Updates after successful step
    if K <= K_old % If step was successful
        K_old = K; % Update cost function to compare with in next iterations
        err_old = err; % Update error to calculate next parameter update
        if nargout >= 3
            try
                model.u_val = []; % Do not save u_val over and over in all successful models
            catch
                % Nothing to do
            end
            model.Cost = K; % Save the obtained cost
            models = [models,model]; % Collect models after successful step
        end
        pack; % Consolidate workspace memory
        model = model_new; % u_val saved in final model (best on estimation data), cost not saved
        states = states_new; %#ok States of the model used in analytical calculation of the Jacobian
        A = model.A; %#ok State matrix used in analytical calculation Jacobian
        C = model.C; %#ok Output matrix used in analytical calculation Jacobian
        E = model.E; %#ok Coefficients state nonlinearity used in analytical calculation Jacobian
        F = model.F; %#ok Coefficients output nonlinearity used in analytical calculation Jacobian
        if sat
            satCoeff = model.satCoeff; %#ok obsolete
        end
    end
end
y_mod = fFilterNLSS(model,u,model.x0,model.u0); % Output of the optimized model

disp('End of L.M. Optimization.')
cd(cwd); % Return to original working directory
end

% ---------------- Help functions ----------------
function [d_powers,d_coeff]=lfPolyDerivative(powers)
%LFPOLYDERIVATIVE Calculate derivative of a multivariate polynomial
d_coeff = permute(powers,[1 3 2]); % Polynomial coefficients of the derivative
n = size(powers,2); % Number of terms
d_powers = repmat(powers,[1 1 n]); % Terms of the derivative
for i = 1:n
    d_powers(:,i,i) = abs(powers(:,i)-1); % Derivative w.r.t. variable i has one degree less in variable i than original polynomial
                                          % If original polynomial is constant w.r.t. variable i, then the derivative is zero, but
                                          % take abs to avoid a power -1 (zero coefficient anyway)
    %{
    % This would be more correct, but is slower
    d_powers(:,i,i) = powers(:,i) - 1;
    d_powers(powers(:,i) == 0,:,i) = 0;
    %}
end
end

function model = lfVec2Par(delta_vector,model)
%LFVEC2PAR Update state-space matrices in model from parameter vector update
% delta_vector is a vector with parameter updates, where the parameters are
% stored as follows: theta = [vec(A.');
%                             vec(B.');
%                             vec(E.');
%                             vec(F.');
%                             vec(satCoeff);
%                             vec(D.');
%                             vec(C.');
%                             x0;
%                             u0]
A = zeros(size(model.A)); % Initialize state matrix update
A(:) = delta_vector(1:model.n^2); % Update in state matrix
model.A = model.A + A'; % Update the state matrix
delta_vector(1:model.n^2) = []; % Remove update in state matrix from update vector

B = zeros(size(model.B))'; % Initialize input matrix update
B(:) = delta_vector(1:model.n*model.m); % Update in input matrix
model.B = model.B + B'; % Update the input matrix
delta_vector(1:model.n*model.m) = []; % Remove update in input matrix from the update vector

E = zeros(size(model.E))'; % Initialize E matrix update
E(model.xactive) = delta_vector(1:length(model.xactive)); % Update in E matrix
model.E = model.E + E'; % Update the E matrix
delta_vector(1:length(model.xactive)) = []; % Remove update in E matrix from the update vector

F = zeros(size(model.F))'; % Initialize F matrix update
F(model.yactive) = delta_vector(1:length(model.yactive)); % Update in F matrix
model.F = model.F + F'; % Update the F matrix
delta_vector(1:length(model.yactive)) = []; % Remove update in F matrix from the update vector

if model.sat
    model.satCoeff = model.satCoeff + delta_vector(1:model.n); % Update saturation coefficients
    delta_vector(1:model.n) = []; % Remove update in saturation coefficients from the update vector
end

D = zeros(size(model.D))'; % Initialize direct feed-through matrix update
D(:) = delta_vector(1:model.p*model.m); % Update in direct feed-through matrix
model.D = model.D + D'; % Update the direct feed-through matrix
delta_vector(1:model.p*model.m) = []; % Remove update in direct feed-through matrix from the update vector

C = zeros(size(model.C))'; % Initialize output matrix update
C(:) = delta_vector(1:model.p*model.n); % Update in output matrix
model.C = model.C + C'; % Update the output matrix
delta_vector(1:model.p*model.n) = []; % Remove update in output matrix from the update vector

x0 = zeros(model.n,1); % Initialize initial state update
x0(model.x0active) = delta_vector(1:length(model.x0active)); % Update in initial state
model.x0 = model.x0 + x0; % Update the initial state
delta_vector(1:length(model.x0active)) = []; % Remove update in initial state from the update vector

u0 = zeros(model.m,1); % Initialize initial state update
u0(model.u0active) = delta_vector(1:length(model.u0active)); % Update in initial state
model.u0 = model.u0 + u0; % Update the initial state
delta_vector(1:length(model.u0active)) = []; % Remove update in initial state from the update vector

% Check compatibility of model and update vector
if ~isempty(delta_vector)
    error('lfVec2Par error: vector is not empty after parsing all parameters')
end
end

%}

%--------------------------------------------------------------------------
% Version 1.0 = Version 1.1
%--------------------------------------------------------------------------
%{

% Set default value for weighting if not specified
if nargin < 5
    W = []; % No weighting
end

% Set default value for Levenberg-Marquardt parameter if not specified
if (nargin < 6) || isempty(lambda)
    lambda = 0; % Later on calculated as the dominant singular value of the Jacobian
end

% Set default value for LambdaJump
% Each LambdaJump iterations, the Levenberg-Marquardt parameter is made
% smaller, so that the algorithm leans more towards a Gauss-Newton
% algorithm (converges faster) and less to a gradient descent method
% (converges in a larger range)
if (nargin < 7) || isempty(LambdaJump)
    LambdaJump = 1001;
end

% Extract parameters for later use from input arguments
N                 = length(u); % Number of samples
[A,C,E,F]         = deal(model.A,model.C,model.E,model.F); %#ok State-space matrices used in analytical calculation Jacobian
[n,m,p]           = deal(model.n,model.m,model.p); %#ok Number of states, inputs, and outputs used in analytical calculation Jacobian
[n_nx,n_ny]       = deal(model.n_nx,model.n_ny); % Number of terms in state and output equation
[xpowers,ypowers] = deal(model.xpowers,model.ypowers); % All possible terms in state and output equation
xactive           = model.xactive; %#ok Active terms in state equation used in analytical calculation Jacobian
yactive           = model.yactive; %#ok Active terms in output equation used in analytical calculation Jacobian
[T1,T2,sat]       = deal(model.T1, model.T2,model.sat); % Transient and saturation parameters

% Extract parameters initial conditions
B = model.B; %#ok Input matrix used in analytical calculation Jacobian
if isfield(model,'x0')
    x0 = model.x0; % Initial state vector
else
    x0 = zeros(n,1); % Zero initial state vector if not specified
    model.x0 = x0;
end
if isfield(model,'u0')
    u0 = model.u0; % Initial input
else
    u0 = zeros(m,1); % Zero initial input if not specified
    model.u0 = u0;
end
if isfield(model,'x0active')
    x0active = model.x0active; %#ok: Indices active elements in initial state vector used in analytical calculation Jacobian
else
    x0active = []; % By default no optimization on initial state vector
    model.x0active = x0active;
end
if isfield(model,'u0active')
    u0active = model.u0active; %#ok: Indices active elements in initial input vector used in analytical calculation Jacobian
else
    u0active = []; % By default no optimization on initial input vector
    model.u0active = u0active;
end

% Compute the (transient-free) modeled output and the corresponding states
[y_mod,states] = fFilterNLSS(model,u,x0,u0); %#ok: states used in analytical Jacobian calculation script

% Determine if weighting is in frequency or time domain: only implemented for periodic signals.
NFD = size(W,3); % Number of frequency bins where weighting is specified (e.g. NFD = floor(Npp/2), where Npp is the number of samples in one period and one phase realization for a multisine excitation)
if isempty(W)
    % No weighting
    FreqWeighting = false; % No frequency weighting
    W = ones(N,p); % No weighting
elseif NFD > 1
    % Frequency weighting
    FreqWeighting = true; % Frequency weighting
    R = round((N-T2)/NFD/2); % Number of realizations
    if mod(N-T2,R) ~= 0
        error('Transient handling and weighting matrix are incompatible')
    end
else
    % Time-domain weighting
    FreqWeighting = false; % No frequency weighting
end

% If T2 is a scalar, it denotes the number of transient points to discard.
% If it is a vector, it denotes the indices of the points to discard, e.g.
% when several data sequences were put together.
without_T2 = fComputeIndicesTransientRemovalArb(T2,N,p); % Samples that are not discarded

% Compute the (weighted) error signal
err_old = y_mod(without_T2) - y(without_T2); % Error signal (without transient for aperiodic signals)
if FreqWeighting
    err_old = reshape(err_old,[(N-T2)/R,R,p]);
    err_old = fft(err_old); % Go to the frequency domain
    err_old = err_old(1:NFD,:,:); % Select only the positive half of the spectrum
    err_old = permute(err_old,[3 1 2]); % p x NFD x R
    err_old = fWeightJacobSubSpace(err_old,W,p,NFD,R); % Add weighting
    K_old = fVec(err_old)'*fVec(err_old); % Calculate cost function
    err_old = permute(err_old,[2 3 1]); % NFD x R x p
    err_old = err_old(:); % NFD R p x 1
    err_old = fReIm(err_old); % Split in real and imaginary part
else
    err_old = err_old.*W(without_T2); % Time-domain weighting
    K_old   = err_old'*err_old; % Compute cost function
end

% Initialization of Levenberg-Marquardt algorithm
Count = 1; % Iteration number
models = []; % Collection of all models after a successful step
Cost = NaN(MaxCount,1); % Sequence of cost functions after a successful step

% Compute the rms value of the weighted output to later on calculate the
% relative error after each successful iteration
if ~FreqWeighting
    rms_y = rms(W(without_T2).*y(without_T2));
end

% Compute the derivatives of the polynomials zeta and eta
[xd_powers,xd_coeff] = lfPolyDerivative(xpowers); %#ok Derivatives of zeta used in analytical calculation Jacobian
[yd_powers,yd_coeff] = lfPolyDerivative(ypowers); %#ok Derivatives of eta used in analytical calculation Jacobian

% Extract the transient part of the input
uT = u(fComputeIndicesTransient(T1,N),:); % Transient part of the input
NT = length(uT); % Length of the transient part

% Should the model be stable on the validation set?
if all(isfield(model,{'u_val','max_out'}))
    stabilisation = true; % Yes, model stable on validation set
else
    stabilisation = false; % No
end

% Prepare for Levenberg-Marquardt optimization
cwd = pwd; % The current working directory
cd(tempdir); % Switch to the system's temporary folder
warning('off','MATLAB:pack:InvalidInvocationLocation'); % Disable warning
pack; % Consolidate workspace memory
disp('Starting L.M. Optimization...')
% inittime = clock; % Save current time to estimate the end time during the optimization 

% Levenberg-Marquardt optimization
while Count < MaxCount
    % Compute the Jacobian w.r.t. the elements in the matrices A, B, E, F, D, and C
    sJacobianAnalytical_x0u0;
    
    % Add frequency- or time-domain weighting
    if FreqWeighting
        numpar = size(J,2); % Number of parameters
        J = reshape(J,[(N-T2)/R,R,p,numpar]); % Npp x R x p x numpar
        J = fft(J); % Go to the frequency domain
        J = J(1:NFD,:,:,:); % Select only the positive half of the spectrum
        J = permute(J,[3 1 2 4]); % p x NFD x R x numpar
        J = reshape(J,[p NFD R*numpar]); % p x NFD x R numpar
        J = fWeightJacobSubSpace(J,W,p,NFD,R*numpar); % Add weighting
        J = reshape(J,[p NFD R numpar]); % p x NFD x R x numpar
        J = permute(J,[2 3 1 4]); % NFD x R x p x numpar
        J = reshape(J,[NFD*R*p,numpar]); % NFD R p x numpar
        J = fReIm(J); % Split in real and imaginary part
    else
        J = J.*repmat(W(without_T2),1,size(J,2)); % Add time-domain weighting
    end

    % (Potentially) improve numerical conditioning
    [J,scaling] = fNormalizeColumns(J); % Normalize columns to unit rms value

    K = K_old; % Initial value of the cost function (=> unsuccessful step)
    [U,S,V] = svd(J,0); % Singular value decomposition of the Jacobian
    clear J;
    pack; % Consolidate workspace memory

    if isnan(K)
        break % Stop nonlinear optimization (typically when initial model is unstable)
    end
    
    nowexit = false; % nowexit is true when a successfull step was found
    
    % Keep looking for a successful step in the parameter update
    while not(nowexit) && Count <= MaxCount % K >= K_old && Count <= MaxCount
        % Estimate the rank of the Jacobian
        s = diag(S); % Singular values of the Jacobian
        tol = NT*eps(max(s)); % Tolerance to determine the number of rank deficiencies
        r = min(sum(s>tol),length(s)-n^2); % Estimate the rank of the Jacobian (there are at least n^2 rank deficiencies)
        
        % Compute parameter update dtheta = -L*err_old, where L = (J.'*J + lambda*eye)^-1*J.' 
        s = diag(s(1:r)./(s(1:r).^2 + lambda^2*ones(r,1))); % Singular values of L
        dtheta = -V(:,1:r)*s*U(:,1:r).'*err_old; % Parameter update
        dtheta = dtheta./scaling'; % Denormalize the estimates (Jacobian was normalized)
        
        % Update the model and save as model_new
        model_new = lfVec2Par(dtheta,model);
        
        % Compute new model output and states
        [y_mod, states_new] = fFilterNLSS(model_new,u,model_new.x0,model_new.u0);
        
        % If the model should be stable on the validation set, ...
        if stabilisation
            % ... compute modeled output on validation set
            model_new.T1 = 0; % No transient handling on validation set
            yStab = fFilterNLSS(model_new,model.u_val,model_new.x0,model_new.u0); % Modeled output
            model_new.T1 = T1; % Transient handling on estimation set
        end
        % If model is not stable on validation set, while this was required, ... 
        if stabilisation && (max(abs(yStab(:))) > model.max_out || any(isnan(yStab(:)))>0)
            K = Inf; % ... set cost function to infinity (=> unsuccessful step), ...
            disp('Unstable for Validation Set');
        else
            % ... else, compute new model error
            err = y_mod(without_T2) - y(without_T2); % Discard transient
            rms_e = rms(err(:)); % Unweighted rms error
            if FreqWeighting
                err = reshape(err,[(N-T2)/R,R,p]); % Npp x R x p
                err = fft(err); % Go to the frequency domain
                err = err(1:NFD,:,:); % Select only the positive half of the spectrum
                err = permute(err,[3 1 2]); % p x Npp x R
                err = fWeightJacobSubSpace(err,W,p,NFD,R); % Add weighting
                K = fVec(err)'*fVec(err); % Weighted least-squares cost function
                err = permute(err,[2 3 1]); % Npp x R x p
                err = err(:); % Npp R p x 1
                err = fReIm(err); % Split in real and imaginary part
            else
                err = err.*W(without_T2); % Time domain weighting (Multiple Output: vec(err))
                K = err'*err; % Cost function
            end
        end
        
        % Check successfullness parameter update
        if K >= K_old % Step was not succesful
            if lambda == 0 % lambda was initialed as zero (lambda = 0 should not be possible later on)
                lambda = S(1,1); % Initial Levenberg-Marquardt parameter = dominant singular value
            else
                lambda = lambda*sqrt(10); % Lean more towards gradient descent method (converges in larger range)
            end;
        elseif isnan(K) % Typically if model is unstable
            lambda = lambda*sqrt(10); % Lean more towards gradient descent method (converges in larger range)
        else % Step was succesful
            lambda = lambda/2; % Lean more towards Gauss-Newton algorithm (converges faster)
            nowexit = true; % Stop searching for successful step
            if FreqWeighting
                % Display weighted mean square cost function (divided by 2*NFD*R*p) and the condition number of the Jacobian 
                disp([num2str(Count) ' Cost Function: ' num2str(K/NFD/R/p/2) ' - Cond: ' num2str(S(1,1)/S(r,r))]);
            else
                % Display rms error normalized with weighted output, the condition number of the Jacobian, the Levenberg Marquardt parameter, and the unweighted rms error
                disp([num2str(Count) '  RelErr ' num2str(rms_e/rms_y) ' - Cond: ' num2str(S(1,1)/S(r,r)) ...
                    ' - Lambda: ' num2str(lambda) '  RmsErr ' num2str(rms_e)]);
            end
            Cost(Count) = rms_e; % Collect unweighted rms errors after each successful step
        end
        
        % Make Levenberg-Marquardt parameter smaller after LambdaJump iterations
        if (rem(Count,LambdaJump)==0 && Count~=0)
            lambda = lambda/10; % Lean more towards Gauss-Newton method (converges faster)
            disp('Forced jump lambda');
        end
        
        % Update the number of iterations
        Count = Count + 1;
    end
    clear U V;
    
    % Updates after successful step
    if K <= K_old % If step was successful
        K_old = K; % Update cost function to compare with in next iterations
        err_old = err; % Update error to calculate next parameter update
        if nargout >= 3
            try
                model.u_val = []; % Do not save u_val over and over in all successful models
            catch
                % Nothing to do
            end
            model.Cost = K; % Save the obtained cost
            models = [models,model]; % Collect models after successful step
        end
        pack; % Consolidate workspace memory
        model = model_new; % u_val saved in final model (best on estimation data), cost not saved
        states = states_new; %#ok States of the model used in analytical calculation of the Jacobian
        A = model.A; %#ok State matrix used in analytical calculation Jacobian
        C = model.C; %#ok Output matrix used in analytical calculation Jacobian
        E = model.E; %#ok Coefficients state nonlinearity used in analytical calculation Jacobian
        F = model.F; %#ok Coefficients output nonlinearity used in analytical calculation Jacobian
        if sat
            satCoeff = model.satCoeff; %#ok obsolete
        end
    end
end
y_mod = fFilterNLSS(model,u,model.x0,model.u0); % Output of the optimized model

disp('End of L.M. Optimization.')
cd(cwd); % Return to original working directory
end

% ---------------- Help functions ----------------
function [d_powers,d_coeff]=lfPolyDerivative(powers)
%LFPOLYDERIVATIVE Calculate derivative of a multivariate polynomial
d_coeff = permute(powers,[1 3 2]); % Polynomial coefficients of the derivative
n = size(powers,2); % Number of terms
d_powers = repmat(powers,[1 1 n]); % Terms of the derivative
for i = 1:n
    d_powers(:,i,i) = abs(powers(:,i)-1); % Derivative w.r.t. variable i has one degree less in variable i than original polynomial
                                          % If original polynomial is constant w.r.t. variable i, then the derivative is zero, but
                                          % take abs to avoid a power -1 (zero coefficient anyway)
    %{
    % This would be more correct, but is slower
    d_powers(:,i,i) = powers(:,i) - 1;
    d_powers(powers(:,i) == 0,:,i) = 0;
    %}
end
end

function model = lfVec2Par(delta_vector,model)
%LFVEC2PAR Update state-space matrices in model from parameter vector update
% delta_vector is a vector with parameter updates, where the parameters are
% stored as follows: theta = [vec(A.');
%                             vec(B.');
%                             vec(E.');
%                             vec(F.');
%                             vec(satCoeff);
%                             vec(D.');
%                             vec(C.');
%                             x0;
%                             u0]
A = zeros(size(model.A)); % Initialize state matrix update
A(:) = delta_vector(1:model.n^2); % Update in state matrix
model.A = model.A + A'; % Update the state matrix
delta_vector(1:model.n^2) = []; % Remove update in state matrix from update vector

B = zeros(size(model.B))'; % Initialize input matrix update
B(:) = delta_vector(1:model.n*model.m); % Update in input matrix
model.B = model.B + B'; % Update the input matrix
delta_vector(1:model.n*model.m) = []; % Remove update in input matrix from the update vector

E = zeros(size(model.E))'; % Initialize E matrix update
E(model.xactive) = delta_vector(1:length(model.xactive)); % Update in E matrix
model.E = model.E + E'; % Update the E matrix
delta_vector(1:length(model.xactive)) = []; % Remove update in E matrix from the update vector

F = zeros(size(model.F))'; % Initialize F matrix update
F(model.yactive) = delta_vector(1:length(model.yactive)); % Update in F matrix
model.F = model.F + F'; % Update the F matrix
delta_vector(1:length(model.yactive)) = []; % Remove update in F matrix from the update vector

if model.sat
    model.satCoeff = model.satCoeff + delta_vector(1:model.n); % Update saturation coefficients
    delta_vector(1:model.n) = []; % Remove update in saturation coefficients from the update vector
end

D = zeros(size(model.D))'; % Initialize direct feed-through matrix update
D(:) = delta_vector(1:model.p*model.m); % Update in direct feed-through matrix
model.D = model.D + D'; % Update the direct feed-through matrix
delta_vector(1:model.p*model.m) = []; % Remove update in direct feed-through matrix from the update vector

C = zeros(size(model.C))'; % Initialize output matrix update
C(:) = delta_vector(1:model.p*model.n); % Update in output matrix
model.C = model.C + C'; % Update the output matrix
delta_vector(1:model.p*model.n) = []; % Remove update in output matrix from the update vector

x0 = zeros(model.n,1); % Initialize initial state update
x0(model.x0active) = delta_vector(1:length(model.x0active)); % Update in initial state
model.x0 = model.x0 + x0; % Update the initial state
delta_vector(1:length(model.x0active)) = []; % Remove update in initial state from the update vector

u0 = zeros(model.m,1); % Initialize initial state update
u0(model.u0active) = delta_vector(1:length(model.u0active)); % Update in initial state
model.u0 = model.u0 + u0; % Update the initial state
delta_vector(1:length(model.u0active)) = []; % Remove update in initial state from the update vector

% Check compatibility of model and update vector
if ~isempty(delta_vector)
    error('lfVec2Par error: vector is not empty after parsing all parameters')
end
end

%}