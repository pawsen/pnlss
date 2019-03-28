function active = fSelectActive(structure,n,m,q,nx)
%FSELECTACTIVE Select active elements in E or F matrix.
%
%	Usage:
%		active = fSelectActive(structure,n,m,q,nx)
%
%	Description:
%		Select the active elements (i.e. those on which optimization will
%		be done) in the E or F matrix. In particular, the linear indices
%		(see also sub2ind and ind2sub) of the active elements in the
%		transpose of the E or F matrix are calculated.
%
%	Output parameters:
%		active : linear indices of the active elements in the transpose of
%		         the E or F matrix
%
%	Input parameters:
%		structure : string indicating which elements in the E or F matrix
%		            are active. The possibilities are 'diagonal',
%		            'inputsonly', 'statesonly', 'nocrossprod', 'affine',
%		            'affinefull', 'full', 'empty', 'nolastinput', or
%		            num2str(row_E). Below is an explanation of each of
%		            these structures:
%                   'diagonal' : active elements in row j of the E matrix
%                                are those corresponding to pure nonlinear
%                                terms in state j (only for state equation)
%                   'inputsonly' : only terms in inputs
%                   'statesonly' : only terms in states
%                   'nocrossprod' : no cross-terms
%                   'affine' : only terms that are linear in one state
%                   'affinefull' : only terms that are linear in one state
%                                  or constant in the states
%                   'full' : all terms
%                   'empty' : no terms
%                   'nolastinput' : no terms in last input (implemented
%                                   since version 1.1)
%                   num2str(row_E) : only row row_E in E matrix is active
%                                    (only for state equation)
%       n : number of states
%       m : number of inputs
%       q : number of rows in corresponding E/F matrix
%           q = n if E matrix is considered,
%           q = p if F matrix is considered
%       nx : degrees of nonlinearity in E/F matrix
%
%	Example:
%		n = 2; % Number of states
%       m = 1; % Number of inputs
%       p = 1; % Number of outputs
%       nx = 2; % Degree(s) of nonlinearity
%       terms = fCombinations(n+m,nx); % Powers of all possible terms in n+m inputs of degree(s) nx
%       % => terms = [2 0 0;
%       %             1 1 0;
%       %             1 0 1;
%       %             0 2 0;
%       %             0 1 1;
%       %             0 0 2];
%       % There are six quadratic terms in the two states x1 and x2, and
%       % the input u, namely x1^2, x1*x2, x1*u, x2^2, x2*u, and u^2.
%       % The matrix E is a 2 x 6 matrix that contains the polynomial
%       % coefficients in each of these 6 terms for both state updates.
%       % The active elements will be calculated as linear indices in the
%       % transpose of E, hence E can be represented as
%       % E = [e1 e2 e3 e4  e5  e6;
%       %      e7 e8 e9 e10 e11 e12];
%       % The matrix F is a 1 x 6 matrix that contains the polynomial
%       % coefficients in each of the 6 terms for the output equation.
%       % The matrix F can be represented as
%       % F = [f1 f2 f3 f4 f5 f6];
%
%       % Diagonal structure
%       activeE = fSelectActive('diagonal',n,m,n,nx);
%       % => activeE = [1 10].';
%       % Only e1 and e10 are active. This corresponds to a term x1^2 in
%       % the first state equation and a term x2^2 in the second state
%       % equation.
%
%       % Inputs only structure
%       activeE = fSelectActive('inputsonly',n,m,n,nx);
%       % => activeE = [6 12].';
%       % Only e6 and e12 are active. This corresponds to a term u^2 in
%       % both state equations. In all other terms, at least one of the
%       % states (possibly raised to a certain power) is a factor.
%       activeF = fSelectActive('inputsonly',n,m,p,nx);
%       % => activeF = 6;
%       % Only f6 is active. This corresponds to a term u^2 in the output
%       % equation.
%
%       % States only structure
%       activeE = fSelectActive('statesonly',n,m,n,nx);
%       % => activeE = [1 2 4 7 8 10].';
%       % Only e1, e2, e4, e7, e8, and e10 are active. This corresponds to
%       % terms x1^2, x1*x2, and x2^2 in both state equations. In all other
%       % terms, the input (possibly raised to a certain power) is a
%       % factor.
%
%       % No cross products structure
%       activeE = fSelectActive('nocrossprod',n,m,n,nx);
%       % => activeE = [1 4 6 7 10 12].';
%       % Only e1, e4, e6, e7, e10, and e12 are active. This corresponds to
%       % terms x1^2, x2^2, and u^2 in both state equations. All other
%       % terms are crossterms where more than one variable is present as a
%       % factor.
%
%       % State affine structure
%       activeE = fSelectActive('affine',n,m,n,nx);
%       % => activeE = [3 5 9 11].';
%       % Only e3, e5, e9, and e11 are active. This corresponds to
%       % terms x1*u and x2*u in both state equations, since in these terms
%       % only one state appears, and it appears linearly. 
%
%       % Full state affine structure
%       activeE = fSelectActive('affinefull',n,m,n,nx);
%       % => activeE = [3 5 6 9 11 12].';
%       % Only e3, e5, e6, e9, e11, and e12 are active. This corresponds to
%       % terms x1*u, x2*u and u^2 in both state equations, since in these
%       % terms at most one state appears, and if it appears, it appears
%       % linearly.
%
%       % Full structure
%       activeE = fSelectActive('full',n,m,n,nx);
%       % => activeE = (1:12).';
%       % All elements in the E matrix are active.
%
%       % Empty structure
%       activeE = fSelectActive('empty',n,m,n,nx);
%       % => activeE = [];
%       % None of the elements in the E matrix are active.
%
%       % One row in E matrix structure
%       row_E = 2; % Select which row in E is active
%       activeE = fSelectActive('row_E',n,m,n,nx);
%       % => activeE = [7 8 9 10 11 12].';
%       % Only the elements in the second row of E are active
%
%       % No terms in last input structure
%       % This is useful in a polynomial nonlinear state-space (PNLSS)
%       % model when considering the initial state as a parameter. The
%       % state at time one can be estimated by adding an extra input
%       % u_art(t) that is equal to one at time zero and zero elsewhere.
%       % Like this, an extended PNLSS model is estimated, where the last
%       % column in its B matrix corresponds to the state at time one in  
%       % the original PNLSS model. To ensure that the optimization is only
%       % carried out on the parameters of the original PNLSS model, only
%       % the corresponding coefficients in the E/F matrix should be
%       % selected as active.
%       terms_extended = fCombinations(n+m+1,nx); % Powers of all possible terms with one extra input
%       % => terms_extended = [2 0 0 0;
%       %                      1 1 0 0;
%       %                      1 0 1 0;
%       %                      1 0 0 1;
%       %                      0 2 0 0;
%       %                      0 1 1 0;
%       %                      0 1 0 1;
%       %                      0 0 2 0;
%       %                      0 0 1 1;
%       %                      0 0 0 2];
%       % The nonlinear terms in the extra input should not be considered
%       % for optimization.
%       activeE_extended = fSelectActive('nolastinput',n,m+1,n,nx);
%       % => activeE_extended = [1 2 3 5 6 8 11 12 13 15 16 18].';
%       % Only the terms where the last input is raised to a power zero are
%       % active. This corresponds to the case where all terms in the
%       % original PNLSS model are active.
%       % The example below illustrates how to combine a certain structure
%       % in the original model (e.g. 'nocrossprod') with the estimation of
%       % the initial state.
%       activeE_extended = fSelectActive('nolastinput',n,m+1,n,nx);
%       activeE_extended = activeE_extended(fSelectActive('nocrossprod',n,m,n,nx));
%       % => activeE_extended = [1 5 8 11 15 18].';
%       % This corresponds to the terms x1^2, x2^2, and u1^2 in both rows
%       % of the E_extended matrix, and thus to all terms in the original
%       % model, except for the crossterms.
%       % Note that an alternative approach is to include the initial state
%       % in the parameter vector (see also fLMnlssWeighted_x0u0).
%
%	Versions:
%		1.0 : November 25, 2015
%       1.1 : December 4, 2015
%           Added structure 'nolastinput'
%       1.2 : April 20, 2016
%           Help updated
%
%	Copyright (c) Vrije Universiteit Brussel – dept. ELEC
%   All rights reserved.
%   Software can be used freely for non-commercial applications only.
%   Disclaimer: This software is provided “as is” without any warranty.
%
%	See also sub2ind, ind2sub, fCombinations, fLMnlssWeighted_x0u0

%--------------------------------------------------------------------------
% Version 1.1 = Version 1.2
%   Added structure 'nolastinput'
%--------------------------------------------------------------------------
% {

combis = fCombinations(n+m,nx); % All possible nonlinear terms of degrees nx in n+m inputs
n_nl = size(combis,1); % Number of terms

switch structure
    case 'diagonal'
        % Diagonal structure requires as many rows in E (or F) matrix as the number of states
        if n ~= q
            error('Diagonal structure can only be used in state equation, not in output equation')
        end
        % Find terms that consist of one state, say x_j, raised to a nonzero power
        active = find(sum(combis(:,1:n) ~= 0,2) == 1 & sum(combis(:,n+1:end) ~= 0,2) == 0);
        % Select these terms only for row j in the E matrix
        for i = 1:length(active) % Loop over these terms
            temp = find(combis(active(i,1),:) ~= 0); % Which state variable is raised to a nonzero power
            active(i) = active(i) + (temp-1)*n_nl; % Linear index of active term in transpose of E
        end
    case 'inputsonly'
        % Find terms where all states are raised to a zero power
        active = find(sum(combis(:,1:n) ~= 0,2) == 0);
        % Select these terms for all rows in E/F matrix
        active = repmat(active,[1,q]) + repmat(0:n_nl:(q-1)*n_nl,[length(active),1]);
        active = active(:);
    case 'statesonly'
        % Find terms where all inputs are raised to a zero power
        active = find(sum(combis(:,n+1:end) ~= 0,2) == 0);
        % Select these terms for all rows in E/F matrix
        active = repmat(active,[1,q]) + repmat(0:n_nl:(q-1)*n_nl,[length(active),1]);
        active = active(:);
    case 'nocrossprod'
        % Find terms where only one variable (state or input) is raised to a nonzero power
        active = find(sum(combis ~= 0,2) == 1);
        % Select these terms for all rows in E/F matrix
        active = repmat(active,[1,q]) + repmat(0:n_nl:(q-1)*n_nl,[length(active),1]);
        active = active(:);
    case 'affine'
        % Find terms where only one state is raised to power one, and all
        % others to power zero. There are no conditions on the powers in
        % the input variables
        active = find(sum(combis(:,1:n) ~= 0,2) == 1 & sum(combis(:,1:n),2) == 1);
        % Select these terms for all rows in E/F matrix
        active = repmat(active,[1,q]) + repmat(0:n_nl:(q-1)*n_nl,[length(active),1]);
        active = active(:);
    case 'affinefull'
        % Find terms where at most one state is raised to power one, and
        % all others to power zero. There are no conditions on the powers
        % in the input variables
        active = find(sum(combis(:,1:n) ~= 0,2) <= 1 & sum(combis(:,1:n),2) <= 1);
        % Select these terms for all rows in E/F matrix
        active = repmat(active,[1,q]) + repmat(0:n_nl:(q-1)*n_nl,[length(active),1]);
        active = active(:);
    case 'full'
        % Select all terms in E/F matrix
        active = (1:q*n_nl).';
    case 'empty'
        % Select no terms
        active = [];
    case 'nolastinput'
        if m > 0
            % Find terms where last input is raised to power zero
            active = find(combis(:,end) == 0);
            % Select these terms for all rows in E/F matrix
            active = repmat(active,[1,q]) + repmat(0:n_nl:(q-1)*n_nl,[length(active),1]);
            active = active(:);
        else
            % If there is no input, then select all terms in E/F matrix
            warning('There is no input')
            active = fSelectActive('full',n,m,q,nx);
        end
    otherwise
        % Check if one row in E is selected
        row_E = str2double(structure);
        if isscalar(row_E) && ismember(row_E,1:n)
            % If so, select that row
            active = (row_E-1)*n_nl + (1:n_nl).';
        else
            error('Unknown structure')
        end
end

% Sort the active elements
active = sort(active);

%}

%--------------------------------------------------------------------------
% Version 1.0
%--------------------------------------------------------------------------
%{

combis = fCombinations(n+m,nx); % All possible nonlinear terms of degrees nx in n+m inputs
n_nl = size(combis,1); % Number of terms

switch structure
    case 'diagonal'
        % Diagonal structure requires as many rows in E (or F) matrix as the number of states
        if n ~= q
            error('Diagonal structure can only be used in state equation, not in output equation')
        end
        % Find terms that consist of one state, say x_j, raised to a nonzero power
        active = find(sum(combis(:,1:n) ~= 0,2) == 1 & sum(combis(:,n+1:end) ~= 0,2) == 0);
        % Select these terms only for row j in the E matrix
        for i = 1:length(active) % Loop over these terms
            temp = find(combis(active(i,1),:) ~= 0); % Which state variable is raised to a nonzero power
            active(i) = active(i) + (temp-1)*n_nl; % Linear index of active term in transpose of E
        end
    case 'inputsonly'
        % Find terms where all states are raised to a zero power
        active = find(sum(combis(:,1:n) ~= 0,2) == 0);
        % Select these terms for all rows in E/F matrix
        active = repmat(active,[1,q]) + repmat(0:n_nl:(q-1)*n_nl,[length(active),1]);
        active = active(:);
    case 'statesonly'
        % Find terms where all inputs are raised to a zero power
        active = find(sum(combis(:,n+1:end) ~= 0,2) == 0);
        % Select these terms for all rows in E/F matrix
        active = repmat(active,[1,q]) + repmat(0:n_nl:(q-1)*n_nl,[length(active),1]);
        active = active(:);
    case 'nocrossprod'
        % Find terms where only one variable (state or input) is raised to a nonzero power
        active = find(sum(combis ~= 0,2) == 1);
        % Select these terms for all rows in E/F matrix
        active = repmat(active,[1,q]) + repmat(0:n_nl:(q-1)*n_nl,[length(active),1]);
        active = active(:);
    case 'affine'
        % Find terms where only one state is raised to power one, and all
        % others to power zero. There are no conditions on the powers in
        % the input variables
        active = find(sum(combis(:,1:n) ~= 0,2) == 1 & sum(combis(:,1:n),2) == 1);
        % Select these terms for all rows in E/F matrix
        active = repmat(active,[1,q]) + repmat(0:n_nl:(q-1)*n_nl,[length(active),1]);
        active = active(:);
    case 'affinefull'
        % Find terms where at most one state is raised to power one, and
        % all others to power zero. There are no conditions on the powers
        % in the input variables
        active = find(sum(combis(:,1:n) ~= 0,2) <= 1 & sum(combis(:,1:n),2) <= 1);
        % Select these terms for all rows in E/F matrix
        active = repmat(active,[1,q]) + repmat(0:n_nl:(q-1)*n_nl,[length(active),1]);
        active = active(:);
    case 'full'
        % Select all terms in E/F matrix
        active = (1:q*n_nl).';
    case 'empty'
        % Select no terms
        active = [];
    otherwise
        % Check if one row in E is selected
        row_E = str2double(structure);
        if isscalar(row_E) && ismember(row_E,1:n)
            % If so, select that row
            active = (row_E-1)*n_nl + (1:n_nl).';
        else
            error('Unknown structure')
        end
end

% Sort the active elements
active = sort(active);

%}