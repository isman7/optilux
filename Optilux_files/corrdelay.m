function [delay,wrn,rho,Iric]=corrdelay(Iric,pat,Nt,Nsymb,opt)

%CORRDELAY System delay by cross-correlation measurement
%   DELAY=CORRDELAY(IRIC,PAT,NT,NSYMB) evaluates the system delay DELAY
%   (in symbols) between the electric signal that drives the transmitter
%   laser and the received signal IRIC. For binary alphabets IRIC is a
%   vector of real, while in multi-level formats IRIC is a complex vector.
%   The transmitted electric signal is artificially created by repeating  
%   the pattern PAT in order to have NT points x symbol (see RESET_ALL). 
%   NSYMB is the number of symbols. 
%
%   E.g. if PAT='1101' and Nt=4 the reference transmitted signal turns out
%   to be non-return to zero equal to:
%   
%   x=1111111100001111
%
%   The cross-correlation is the convolution between the vector conj(x) and 
%   IRIC.
%   
%   This function measures the delay as the coordinate of the maximum of
%   the cross-correlation. If the  difference between the largest and the 
%   second maximum is lower than a threshold (initialized at the beginning 
%   of this function), the function display a warning of low accuracy. 
%   Low accuracy is an indicator that the delay may be wrong. 
%
%   DELAY=CORRDELAY(IRIC,PAT,NT,NSYMB,OPT) has the optional flag OPT.
%   If OPT='phase' CORRDELAY treats IRIC as a complex signal (e.g. QPSK 
%   signal). 
%
%   [DELAY,WRN,RHO]=CORRDELAY(IRIC,PAT,NT,NSYMB,OPT) also returns the 
%   correlation coefficient RHO of the signal and reference one. 
%
%   WRN is a flag equal to true if the the delay is measured with
%   low accuracy.
%   Low accuracy may happen in presence of amplified spontaneous emission 
%   noise or with big distortions. 
%   In such cases, if possible, it is better to use the theoretical delay 
%   saved into global variable GSTATE.DELAY.
%
%   [DELAY,WRN,RHO,IRICA]=CORRDELAY(IRIC,PAT,NT,NSYMB,OPT) returns in IRICA 
%	the angle of the input IRIC after removing phase shift ambiguity (useful 
%	in multi-level formats).
%
%   See Also EVAL_EYE, DSP4COHDEC
%
%   Author: Paolo Serena, 2008
%   University of Parma, Italy

%    This file is part of Optilux, the optical simulator toolbox.
%    Copyright (C) 2008  Paolo Serena, <serena@tlc.unipr.it>
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

if nargin < 5
    isopt = false;
elseif strcmp(opt,'phase')
    isopt = true;
else
    error('wrong flag for opt.');
end

MINERR = 1e-4; % [dB].relerr < MINERR is an indicator of a periodic pattern
MAXERR = 0.5;  % [dB].relerr < MAXERR is an indicator that the delay is wrong
NPHI = 71;    % number of phases under for removing phase ambiguity

wrn=false;
Nfft = Nsymb*Nt;
if isopt
    allphi = (0:NPHI-1)*2*pi/(floor(NPHI));% linspace except 2*pi
    refsig = reshape(repmat(fastexp(pat),1,Nt).',Nfft,1);
    crosscorr = (ifft( fft(Iric) .* conj(fft(refsig))));
    vettmax = zeros(1,NPHI);
    vettpos = zeros(1,NPHI);
    for kk=1:NPHI
    	temp = cos(allphi(kk))*real(crosscorr) - sin(allphi(kk))*imag(crosscorr);
    	[vettmax(kk),vettpos(kk)]=max(temp);
    end
    [maxc,posmax] = max(vettmax); 
    posc = vettpos(posmax);
    Iric = angle(Iric*fastexp(-allphi(posmax))); % remove phase for eye evaluation
    crosscorr = real(crosscorr*fastexp(allphi(posmax)));
else
    refsig = reshape(repmat(pat,1,Nt)',Nfft,1);
    % refsig is an ideal reference signal for delay evaluation
    crosscorr = real(ifft( fft(Iric) .* conj(fft(refsig))));
    [maxc,posc] = max(crosscorr);
end

delay = nmod(posc,Nfft)-1; % system delay

% Now check if the delay is correct

ii=find(diff(sign(diff(crosscorr))) == -2) + 1;% all the maxima of crosscorr 
[allmax,posall]=sort(crosscorr(ii),'descend'); 
if ii(posall(1)) == posc
    inderr = 2; % use the second largest maximum
else
    inderr = 1; % means that the maximum is at one boundary of crosscorr
end
relerr = 10*log10(crosscorr(posc)/allmax(inderr));  % relative error
if (relerr > MINERR) && (relerr < MAXERR)
    warning('optilux:corrdelay','low accuracy in delay evaluation');
    wrn=true;
end
delay = (delay+Nt/2)/Nt;    % delay in bits. Double number.
rho = maxc/Nfft*2;
% Note: +Nt/2 because the first bit is centered on the first sample.
