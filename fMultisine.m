function [u,lines,non_exc_odd_or_all,non_exc_even] = fMultisine(N,kindOrLines,M,R,Nblock)
%FMULTISINE Generate a random-phase multisine signal with a flat spectrum.
%
%	Usage:
%		[u,lines,non_exc_odd,non_exc_even] = fMultisine(N,kind,M,R,Nblock)
%		[u,lines,non_exc_odd,non_exc_even] = fMultisine(N,lines,M,R,Nblock)
%		[u,lines,non_exc] = fMultisine(N,kind,M,R,Nblock)
%		[u,lines,non_exc] = fMultisine(N,lines,M,R,Nblock)
%       u = fMultisine()
%
%	Description:
%		fMultisine generates R realizations of a random-phase multisine
%		signal with N samples and excitation up to line M (line 1, i.e. DC,
%		not included). The excited lines in this band can be specified
%		either as a vector with the excited lines or as a string indicating
%		the type of multisine (e.g. all lines excited, only odd lines
%		excited, etc. (see info with input parameter kindOrLines)).
%
%	Output parameters:
%       u : N x R matrix with multisine signal(s)
%       lines : Excited lines
%       non_exc_odd_or_all : Non-excited odd lines (if there are four
%                            output arguments) or all non-excited lines (if
%                            there are three output arguments)
%       non_exc_even : Non-excited even lines
%
%	Input parameters:
%		N : Number of samples (optional, default = 1024)
%       kindOrLines : String indicating the type of multisine: 
%                     'full': all lines excited
%                     'odd': only odd lines excited
%                     'special odd': odd frequencies 8*k+1 and 8*k+3 (=>
%                                    lines 8*k+2 and 8*k+4) excited for 
%                                    k = 0, 1, ...
%                     'random odd': odd lines excited except for one random
%                                   detection line in each group of Nblock
%                                   consecutive odd lines
%                     or a vector with the excited lines
%                     (optional, default = 'full')
%       M : Last excited line (optional, default = 90% of Nyquist)
%       R : Number of realizations (optional, default = 1)
%       Nblock : Block length for random odd multisine (optional, default = 4)
%
%	Examples:
%       % Generate two realizations of a full multisine with 1000 samples
%       % and excitation up to one third of the Nyquist frequency
%       N = 1000; % One thousand samples
%       kind = 'full'; % Full multisine
%       M = floor(N/6); % Excitation up to one sixth of the sample frequency
%       R = 2; % Two phase realizations
%       u = fMultisine(N,kind,M,R); % Multsine signals
%       % Check spectra
%       figure
%       subplot(2,2,1)
%       plot(db(fft(u(:,1))),'+')
%       xlabel('Frequency line')
%       ylabel('Amplitude (dB)')
%       title({'Phase realization 1'})
%       subplot(2,2,2)
%       plot(db(fft(u(:,2))),'+')
%       xlabel('Frequency line')
%       ylabel('Amplitude (dB)')
%       title({'Phase realization 2'})
%       subplot(2,2,3)
%       plot(angle(fft(u(:,1))),'+')
%       xlabel('Frequency line')
%       ylabel('Phase (rad)')
%       title({'Phase realization 1'})
%       subplot(2,2,4)
%       plot(angle(fft(u(:,2))),'+')
%       xlabel('Frequency line')
%       ylabel('Phase (rad)')
%       title({'Phase realization 2'})
%       % The two realizations have the same amplitude spectrum, but
%       % different phase realizations (uniformly distributed between -pi
%       % and pi)
%
%       % Generate a random odd multisine where the excited odd lines are
%       % split in groups of three consecutive lines and where one line is
%       % randomly chosen in each group to be a detection line (i.e.
%       % without excitation)
%       N = 1000; % One thousand samples
%       kind = 'random odd'; % Random odd multisine
%       M = floor(N/6); % Excitation up to one sixth of the sample frequency
%       R = 1; % One phase realization
%       Nblock = 3; % One out of three consecutive odd lines is randomly selected to be a detection line
%       [u1,lines] = fMultisine(N,kind,M,R,Nblock); % Multisine signal and excited lines
%       % Generate another phase realization of this multisine with the
%       % same excited lines and detection lines
%       u2 = fMultisine(N,lines,[],1); % One realization of a multisine with the same excited lines as u1
%
%       % Change the coloring and rms level of the generated multisine
%       u = fMultisine(); % Default multisine
%       [b,a] = cheby1(2,10,2*0.1); % Filter coefficients
%       U = fft(u); % Multisine spectrum
%       U_colored = freqz(b,a,2*pi*(0:length(u)-1).'/length(u)).*U; % Filtered spectrum
%       u_colored = real(ifft(U_colored)); % Colored multisine signal
%       u_scaled = 2*u_colored/std(u_colored); % Scaled signal to rms value 2 (u_colored is zero-mean)
%
%	Versions:
%		1.0 : March 7, 2016
%       1.1 : April 20, 2016
%           Help updated
%
%	Copyright (c) Vrije Universiteit Brussel – dept. ELEC
%   All rights reserved.
%   Software can be used freely for non-commercial applications only.
%   Disclaimer: This software is provided “as is” without any warranty.

%--------------------------------------------------------------------------
% Version 1.0 = Version 1.1
%--------------------------------------------------------------------------
% {

% Set default values if not specified
if (nargin < 1) || isempty(N)
    N = 1024; % Number of samples
end
if (nargin < 2) || isempty(kindOrLines)
    kind = 'full'; % Kind of multisine
else
    if isa(kindOrLines,'char') % Kind of multisine specified
        kind = kindOrLines;
    else % Lines specified
        kind = 'dummy';
        lines = kindOrLines;
    end
end
if (nargin < 3) || isempty(M)
    M = floor(0.9*N/2); % Last excited line
end
if (nargin < 4) || isempty(R)
    R = 1; % Number of realizations
end
if (nargin < 5) || isempty(Nblock)
    Nblock = 4; % Block length for a random-odd multisine
end

% Set the excited lines
switch lower(kind)
    case 'full'
        lines = 2:M;
    case 'odd'
        lines = 2:2:M;
    case {'specialodd', 'special odd','special'}
        lines = sort([2:8:M,4:8:M]);
    case {'randomodd','random odd','random'}
        lines = 2:2:M;
        Nblocks = floor(length(lines)/Nblock);
        indices = 1:Nblock:Nblocks*Nblock; % Start indices of blocks ...
        indices = indices+floor(Nblock*rand(1,Nblocks)); % ... plus random number
        lines(indices) = []; % Eliminate selected detection lines
    case 'dummy'
        % Do nothing
    otherwise
        error('Unknown type');
end

% Generate multisine signal
U = zeros(N,R);
U(lines,:) = exp(2*pi*1j*rand(length(lines),R));
u = 2*real(ifft(U));

% Non-excited lines
if nargout > 2
    non_lines = 1:N/2;
    non_lines(lines) = [];
    non_odd = non_lines(logical(mod(non_lines-1,2)));
    non_even = non_lines(~mod(non_lines-1,2));
    if nargout > 3
        non_exc_odd_or_all = non_odd;
        non_exc_even = non_even;
    else
        non_exc_odd_or_all = non_lines;
    end
end
%}