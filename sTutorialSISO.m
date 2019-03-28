%STUTORIALSISO Tutorial that shows how to simulate and estimate a SISO PNLSS model.
%
%   Description:
%       sTutorialSISO contains a tutorial that shows the workflow of
%       modeling a single input single output (SISO) polynomial nonlinear
%       state-space (PNLSS) model.
%       The tutorial shows how to 
%           1) generate a multisine signal
%           2) compute the response of a SISO PNLSS model to the multisine
%              input
%           3) estimate a nonparametric linear model from the input and 
%              noisy output data
%           4) estimate a parametric linear state-space model on the
%              nonparametric model
%           5) estimate the parameters in the full PNLSS model
%
%	Copyright (c) Vrije Universiteit Brussel – dept. ELEC
%   All rights reserved.
%   Software can be used freely for non-commercial applications only.
%   Disclaimer: This software is provided “as is” without any warranty.

%% Clear workspace, command window, and figures

clear
clc
close all

%% Add folders to search path

addpath(pwd) % Add the current folder to the search path

%% Load a PNLSS model (will be used to generate data for this example code)

load dTutorialSISO.mat model
true_model = model;
clear model % Will be estimated

%% Generate multisine input (u)

RMSu = 0.05; % Root mean square value for the input signal
N = 1024;    % Number of samples
R = 4;       % Number of phase realizations (one for validation and one for performance testing)
P = 3;       % Number of periods
kind = 'Odd';           % 'Full','Odd','SpecialOdd', or 'RandomOdd': kind of multisine
M = round(0.9*N/2);     % Last excited line
[u,lines,non_exc_odd,non_exc_even] = fMultisine(N, kind, M, R); % Multisine signal, excited and detection lines
u = u/rms(u(:,1))*RMSu; % Scale multisine to the correct rms level
u = repmat(u,[1 1 P]);  % N x R x P
u = permute(u,[1 3 2]); % N x P x R

%% Generate output (y)

% Calculate output of PNLSS system
true_model.T1 = [N 1+(0:P*N:(R-1)*P*N)]; % To generate steady state data
y = fFilterNLSS(true_model,u(:));    % Supposed to be the true system output (simulation)

% Add colored noise to the output
y = reshape(y,[N P R]); % N x P x R
noise = 1e-3*std(y(:,end,end))*randn(size(y)); % Output noise signal
noise(1:end-1,:,:) = noise(1:end-1,:,:) + noise(2:end,:,:); % Do some filtering
y = y + noise; % Noise added to the output

%% Separate the data in estimation, validation, and test set

% Last realization, last period for performance testing
utest = u(:,end,R); utest = utest(:);
ytest = y(:,end,R); ytest = ytest(:);

% One but last realization, last period for validation and model selection
uval = u(:,end,R-1); uval = uval(:);
yval = y(:,end,R-1); yval = yval(:);

% All other realizations for estimation
R = R-2;
u = u(:,:,1:R);
y = y(:,:,1:R);

%% Estimate nonparametric linear model (BLA)

u = permute(u,[1,4,3,2]); % N x m x R x P
y = permute(y,[1,4,3,2]); % N x p x R x P
covY = fCovarY(y); % Noise covariance (frequency domain)

U = fft(u); U = U(lines,:,:,:); % Input spectrum at excited lines
Y = fft(y); Y = Y(lines,:,:,:); % Output spectrum at excited lines

% Estimate best linear approximation, total distortion, and noise distortion
[G,covGML,covGn] = fCovarFrf(U,Y); % G: FRF; covGML: noise + NL; covGn: noise (all only on excited lines)

%% Estimate linear state-space model (frequency domain subspace)

% Choose model order
na   = 2;  % Model order
maxr = 10; % Subspace dimensioning parameter
freq = (lines-1)/N; % Excited frequencies (normalized)

% covGML = repmat(eye(1),[1 1 length(lines)]); % Uncomment for uniform weighting (= no weighting)
models = fLoopSubSpace(freq,G,covGML,na,maxr,100); % All estimated subspace models

% Extract linear state-space matrices from best model on validation data
Nval = length(uval); % Number of samples in validation data
fs = 1; % Sampling frequency
tval = (0:Nval-1)/fs; % Time vector validation data
min_err = Inf; % Initialize minimum error
min_na = NaN; % Initialize model order of best model
for n = na % Loop over model order(s)
    model = models{n}; % Select subspace model of the correct order
    A = model{1}; B = model{2}; C = model{3}; D = model{4}; % Extract state-space matrices
    [A,B,C] = dbalreal(A,B,C); % Compute balanced realization
    yval_hat = lsim(ss(A,B,C,D,1/fs),uval,tval); % Modeled output
    err = yval - yval_hat; % Error signal
    err = sqrt(mean(err(end-N+1:end).^2)); % Rms value of the last period of the error signal
    if err < min_err % If the model is the best so far
        min_na = n; % Update model order of the best model
        min_err = err; % Update minimum error
    end
end
model = models{min_na}; % Select the best model
[A,B,C,D] = model{:}; % Extract all the matrices of the best model
[A,B,C] = dbalreal(A,B,C); % Balanced realization

%% Estimate PNLSS model

% Estimation data
u = mean(u,4); % Average over periods
y = mean(y,4); % Average over periods (be careful that the data are truly steady state)
m = 1; % Number of inputs
p = 1; % Number of outputs
u = u(:); % Concatenate the data: N*P*R x m
y = y(:); % Concatenate the data: N*P*R x p

% Transient settings
NTrans = N; % Add one period before the start of each realization
T1 = [NTrans 1+(0:N:(R-1)*N)]; % Number of transient samples and starting indices of each realization
T2 = 0; % No non-periodic transient handling

% Nonlinear terms
nx = [2 3]; % Nonlinear degrees in state update equation
ny = [2 3]; % Nonlinear degrees in output equation
whichtermsx = 'full'; % Consider all monomials in the state update equation
whichtermsy = 'full'; % Consider all monomials in the output equation

% Settings Levenberg-Marquardt optimization
MaxCount = 100; % Number of Levenberg-Marquardt optimizations
lambda = 100; % Starting value Levenberg-Marquardt damping factor

% Choose model order
n = min_na;

% Initial linear model in PNLSS form
model = fCreateNLSSmodel(A,B,C,D,nx,ny,T1,T2); % Initialize PNLSS model with linear state-space matrices

% Set which monomials will be optimized
model.xactive = fSelectActive(whichtermsx,n,m,n,nx); % Select active monomials in state equation
model.yactive = fSelectActive(whichtermsy,n,m,p,ny); % Select active monomials in output equation

% Protect for unstable models
model.uval = uval; % Optionally, a validation input signal can be passed for which the output of the estimated model should remain bounded
model.max_out = 1000*max(abs(yval)); % Bound on the output of the estimated model (if this bound is violated, the parameter update is considered as unsuccessful)

% Compute linear model error
y_lin = fFilterNLSS(model,u); % Output of the initial linear model on the estimation data
modellintest = model; modellintest.T1 = NTrans; % Change transient parameter for linear model on test data
ytest_lin = fFilterNLSS(modellintest,utest); % Output of the initial linear model on the test data
err_lin = ytest - ytest_lin; % Output error initial linear model on test data

% Set weighting
for kk = size(covY,3):-1:1
    W(:,:,kk) = fSqrtInverse(covY(:,:,kk)); % Frequency weighting
end;
%     W = []; % Uncomment for uniform (= no) weighting

% Levenberg-Marquardt optimization
[model, y_mod, models] = fLMnlssWeighted(u,y,model,MaxCount,W,lambda); % The actual optimisation of the PNLSS model

% Time-domain plot on the estimation data
figure; plot([y y-y_lin y-y_mod]) % Estimation results: output, linear error and PNLSS error
xlabel('Time index')
ylabel('Output (errors)')
legend('output','linear error','PNLSS error')
title('Estimation results')
disp(' ')
disp(['rms(y-y_mod) = ' num2str(rms(y-y_mod)) ' (= Output error of the best PNLSS model on the estimation data)'])
disp(['rms(noise(:))/sqrt(P) = ' num2str(rms(noise(:))/sqrt(P)) ' (= Noise level)'])

% Search best model over the optimisation path on a fresh set of data
valerrs = [];
for i = 1:length(models)
    models(i).T1 = NTrans; % Only one realization in the validation data        
    yval_mod = fFilterNLSS(models(i),uval);  % Output model i on validation data
    valerr = yval - yval_mod; % Output error model i on validation data
    valerrs = [valerrs; rms(valerr)]; % Collect output errors of all models in a vector
end
figure;
plot(db(valerrs));
xlabel('Successful iteration number')
ylabel('Validation error [dB]')
title('Selection of the best model on a separate data set')
[min_valerr,i] = min(valerrs); % Select the best model on the validation data to avoid overfitting
hold on
plot(i,db(min_valerr),'r.','Markersize',10)
model = models(i);

%% Result on test data

% Compute output error on test data
valerr = ytest - fFilterNLSS(model,utest);

% Frequency-domain plot on the test data
figure;
fs = 1; % Just a normalised discrete-time simulation
freq = (0:N-1)/N*fs; % Normalized frequency vector
plottime = [ytest err_lin valerr]; % Test output and output errors linear and PNLSS model (in time domain)
plotfreq = fft(plottime); % Test output and output errors linear and PNLSS model (in frequency domain)
plot(freq(1:end/2),db(plotfreq(1:end/2,:)),'.')
xlabel('Frequency (normalised)')
ylabel('Output (errors) [dB]')
legend('Output','Linear error','PNLSS error')
title('Test results')