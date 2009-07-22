function Eout=rz_carver( Ein, style, options)
%RZ_CARVER carves the optical field to produce return to zero pulses
%   E=RZ_CARVER(E,STYLE, OPTIONS) carves the optical field E with a  
%   Mach-Zehnder interferometer to produce return to zero (RZ) pulses [1]. 
%   The parameter STYLE is a string that designates the type of carving 
%   scheme. Supported schemes are:
%
%       - RZ 33% (STYLE = 'RZ33')
%       - RZ 50% (STYLE = 'RZ50')
%       - RZ 66% (CSRZ, [1]) (STYLE = 'RZ66')
%
%   OPTIONS is an optional structure whose fields can be:
%
%       - exratio: extinction ratio [dB] (default = inf)
%       - bias: bias of the carver (default = 0)
%       - amplitude: Vpi of the carver (default = 1)
%       - nochirp: reduce effect of chirp due to finite exratio [2]
%         (default = false)
%   
%   See also MZ_MODULATOR, LINEAR_MODULATOR, PHASE_MODULATOR, QI_MODULATOR
%   
%   [1] P. J. Winzer, R. J. Essiambre, "Advanced modulation formats for
%   high-capacity optical transport networks," Journal Lightw. Technol.,
%   vol. 24, pp.4711-4728, Dec. 2006.
%   [2] K. Hoon and A. H. Gnauck, "Chirp characteristics of dual-drive 
%   Mach-Zehnder modulator with a finite DC extinction ratio," IEEE Photon. 
%   Technol. Lett., vol. 14, pp. 298-300, Mar. 2002.
%
%   Author: Marco Bertolini, 2009
%   University of Parma, Italy

%    This file is part of Optilux, the optical simulator toolbox.
%    Copyright (C) 2009  Paolo Serena, <serena@tlc.unipr.it>
%			 
%    Optilux is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.
%
%    Optilux is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

global GSTATE;

bias      = 0;
amplitude = 1;
exratio     = inf;
nochirp   = false;

switch style
    case 'RZ50'
        modsig = elecsrc(GSTATE.NT/2,reshape([1 0]'*ones(1,GSTATE.NSYMB),GSTATE.NSYMB*2,1),'cosroll',1,1,'');
    case 'RZ33'
        modsig = fastshift(elecsrc(GSTATE.NT,reshape([1 -1]'*ones(1,GSTATE.NSYMB/2),GSTATE.NSYMB,1),'cosroll',1,1,''),GSTATE.NT/2);
        bias = 1;
    case 'RZ66'
        modsig = elecsrc(GSTATE.NT,reshape([1 -1]'*ones(1,GSTATE.NSYMB/2),GSTATE.NSYMB,1),'cosroll',1,1,'');
    otherwise
        error('error in rz_carver: unsupported carving scheme');
end

if exist('options','var')
    checkfields(options,{'bias','amplitude','exratio','nochirp'});
    if isfield(options,'bias');
        bias = options.bias;
    end
    if isfield(options,'amplitude');
        amplitude = options.amplitude;
    end
    if isfield(options,'exratio');
        exratio = options.exratio;
    end  
    if isfield(options,'nochirp');
        nochirp = options.nochirp;
    end      
end

exr_lin = 10^(-exratio/10);
gamma = (1-sqrt(exr_lin))/(sqrt(exr_lin)+1);

if nochirp
    Phi_U = pi/(1+gamma^2) * (modsig * amplitude + bias);
    Phi_L = -gamma^2*pi/(1+gamma^2) * (modsig * amplitude + bias);    
else
    Phi_U = pi/2 * (modsig * amplitude + bias);
    Phi_L = -pi/2 * (modsig * amplitude + bias);    
end


% Any phase shift due to the interferometric structure is neglected for the
% sake of simplicity
Eout  = j*Ein.*(fastexp(Phi_L)-gamma*fastexp(Phi_U))/(1+gamma);

function elec = elecsrc(Nsps,pat,ptype,duty,roll,inpow,ord)

if max(duty) > 1
    error('error in rz_carver.m: duty must be <= 1');
end
if exist('roll','var') ~= 1
    error('error in rz_carver.m: Undefined roll');
end
if strcmp(ptype,'cosroll') && (max(roll) > 1)
    error('error in interz_carvernsitymod.m: roll-off must be <= 1');
end

if ~(strcmp(ptype,'cosroll') || strcmp(ptype,'dirac') || strcmp(ptype,'idpulse')...
        || strcmp(ptype,'sech') || strcmp(ptype,'tanh'))
    
     flag = 1;   % after modulation, apply the filter.
     ptypet = 'idpulse';     % temp ptype
else
     flag = 0;   % after modulation, exit the function.
     ptypet = ptype;
end


% Modulate
Nsymb = length(pat);
Nfft = Nsymb*Nsps;

elec = zeros(Nfft,1);

elpulse = pulseshape(ptypet,roll,duty,Nsps);  % single pulse

nstart = Nfft-Nsps+1;      % The first elpulse starts at the end of
% the sequence (cyclic periodicity).
nend = Nfft;
elec(nstart:nend) = pat(1)*elpulse(1:Nsps);
elec(1:Nsps) = pat(1)*elpulse(Nsps+1:Nsps*2);
for kbit=2:Nsymb          % all other bits: 2 -> Nsymb
    nstart = (kbit-2)*Nsps+1;
    nend = kbit*Nsps;
    elec(nstart:nend) = elec(nstart:nend)+pat(kbit)*elpulse; % add the kbit-pulse
end

if flag == 1    % filter the signal
    Hf = myfilter(ptype,fftshift(-Nsps/2:1/Nsymb:Nsps/2-1/Nsymb).',roll,ord);
    elec = ifft(fft(elec).* Hf);
    delay = evaldelay(ptype,roll);
    elec=circshift(elec,round(-delay*Nsps));
end

if strcmp(inpow,'power')
    elec = sqrt(abs(elec));
end

%------------------------------------------------
function y=pulseshape(ptype,roll,duty,Nt)

%PULSESHAPE Creates the fundamental pulse
%   Y=PULSESHAPE(PTYPE,ROLL,DUTY,NT) returns in Y a vector [2*NT,1]
%   containing the fundamental pulse whose type is defined in ptype (see
%   intensitymod.m for available types). NT is the number of points x bit,
%   DUTY the bit duty-cycle. ROLL is the pulse roll-off if used, otherwise 
%   any number. The length of Y is 2*NT because the roll-off spreads the 
%   pulse outside the bit time.

elpulse = zeros(Nt*2,1);     % elementary pulse (over two bit times
                             % because the roll-off spreads the pulse).
switch ptype
    case 'cosroll'

        nl = round(0.5*(1-roll)*duty*Nt);    % start index of cos roll-off
        nr = duty*Nt-nl-1;                   % end index of cos roll-off

        nmark = 1:nl;                       % indeces where the pulse is 1
        ncos  = nl:nr;                      % transition region of cos roll-off

        elpulse(Nt+nmark) = 1;
        hperiod = duty*Nt-2*nl;
        if hperiod ~= 0
            elpulse(ncos+Nt+1) = 0.5*(1+cos(pi/(hperiod)*(ncos-nl+0.5)));
        end
        elpulse(1:Nt) = flipud(elpulse(Nt+1:Nt*2)); % first half of the pulse

    case 'dirac'
        elpulse(Nt+1) = 1;
    case 'idpulse'
        nl = round(0.5*duty*Nt);
        if nl == 0
            elpulse(Nt+1) = 1;  % same as Dirac's delta. Why this? Because in this
                                % way you can create a gaussian pulse, for
                                % instance, by filtering this train of delta.
        else
            elpulse(Nt-nl+1:Nt+nl) = 1;
        end
    otherwise
        error('error in rz_carver.m: the pulse ptype does not exist');
end

y=elpulse;
