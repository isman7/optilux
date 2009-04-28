function wrn=reset_all(Nsymb,Nt,Nch,opt1,opt2)

%RESET_ALL Reset all global variables and initialize the simulation.
%   THIS FUNCTION MUST BE CALLED IN EACH SIMULATION.
% 
%   RESET_ALL(NSYMB,NT,NCH) initializes the global struct variable GSTATE.
%   The fields of the structure GSTATE are the following:  
%       
%       GSTATE.NSYMB:   Number of symbols = NSYMB. For binary transmissions
%                       NSYMB is the number of bits.
%       GSTATE.NT:      Number of discrete points x symbol = NT
%       GSTATE.NCH:     Number of channels = NCH
%       GSTATE.FN:      Frequency normalized to the symbol rate (R),
%                       i.e. GSTATE.FN = freq/R with freq the frequency
%                       in [Hertz], R the symbol rate in [baud]. 
%                       
%                                   1/R*DUTY
%                                   <------->
%                                    -------
%                                   |       |
%                                   |       |   single symbol
%                                   |       |
%                                   |       |
%                             -------       -------     
%                             <------------------->
%                               symbol time = 1/R
%
%                       * The LOWEST discrete frequency (resolution) is
%                         1/NSYMB.
%                       * The LARGEST discrete frequency (Nyquist frequency)
%                         is  NT/2.
%                       * The frequencies into GSTATE.FN:
%                             fftshift(-NT/2:1/NSYMB:NT/2-1/NSYMB).
%                       * GSTATE.FN is dimensionless.
%
%       GSTATE.SYMBOLRATE: Symbol rate (R) in [GBaud] equal for all channel
%                       See the examples for how to manage different
%                       channel symbol rates. The symbolrate will be
%                       initialized to a numerical value in ELECTRICSOURCE.
%       GSTATE.PRINT:   true/false. True: Print functions details to file.
%
%
%   RESET_ALL(NSYMB,NT,NCH,OPT1) creates the output directory OPT1 that
%   will collect a summary of each function operation within the file
%   simul_out. The file simul_out is appended each call. The output 
%   directory will also contain any signal printed to file by optilux.  
%   In presence of OPT1 there is the additional field:
%
%       GSTATE.DIR = OPT1
%
%   WRN=RESET_ALL(NSYMB,NT,NCH,OPT1) sets WRN=true if the dimension of
%   simul_out exceeds 50 Mbytes, otherwise WRN=false.
%
%   RESET_ALL(NSYMB,NT,NCH,OPT1,OPT2) with OPT2='noprint' creates the
%   output directory OPT1 but all functions called by optilux will not 
%   print any detail in simul_out.
%
%   RESET_ALL initializes the fundamental constants (like Planck's, speed
%   of light, etc) in the global variable CONSTANTS.
%
%   There are other fields of GSTATE that are set to empty by RESET_ALL and
%   will be initialized by CREATE_FIELD. They are:
%
%       GSTATE.FIELDX: x-component of the electric field. See CREATE_FIELD.
%       GSTATE.FIELDY: y-component of the electric field. See CREATE_FIELD. 
%       GSTATE.FIELDX_TX: copy of GSTATE.FIELDX generated by CREATE_FIELD.  
%                         Useful for back-to-back measurements.
%       GSTATE.FIELDY_TX: same as GSTATE.FIELDX_TX, but for y polarization.
%
%       GSTATE.DELAY:  Overall delay cumulated in the optical line, normalized
%                      to 1/R. Size: [2,NCH] if the y-component is empty, else 
%		       [1,NCH].
%       GSTATE.DISP:    cumulated dispersion [ps/nm] in the system.
%                       Size: same as GSTATE.DELAY.
%
%   Other fields of GSTATE will be initialized by LASERSOURCE. They are:
%
%       GSTATE.LAMBDA:  Channels wavelength [nm]. Size: [1,NCH]
%       GSTATE.POWER:   Transmitted Signal PEAK power [mW]. Size: [1,NCH]
%
%
%   See also: CREATE_FIELD, ELECTRICSOURCE, LASERSOURCE
%
%   Author: Paolo Serena, 2009
%   Contributed by: Armando Vannucci, 2009
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

global CONSTANTS;  % CONSTANTS is a global structure variable.
CONSTANTS.CLIGHT= 299792458;        % speed of light in vacuum [m/s]
CONSTANTS.HPLANCK= 6.62606896e-34;  % Planck's constant [J*s] (CODATA value, 
                                    % std. uncertainty= 3.3e-43)
                                    % NOTE: in October 2005, the National 
                                    % Physical Laboratory reports 6.62607095e-34
CONSTANTS.ECHARGE= 1.602176487e-19;  % Electron's charge [C] (CODATA value, year 2006)
CONSTANTS.KBOLTZMANN= 1.3806504E-23; % Boltzmann's constant [J/oK] (CODATA value, 
                                     % std. uncertainty= 2.4e-31)

global GSTATE;  % GSTATE is a global structure variable.

%%%%%%%%%%%%%%%
MAXBYTES = 50e6;    % warning dimension for simul_out

error(nargchk(3,5,nargin));

switch nargin
    case 3
        GSTATE.PRINT = false; % print on output file simul_out
    case 4
        if ~ischar(opt1)
            error('directory name must be a string');
        end
        if strcmp(opt1,'noprint')
            error('The output directory cannot be called ''noprint''');
        end
        GSTATE.DIR = opt1;
        GSTATE.PRINT = true;
    case 5
        if strcmp(opt1,'noprint')
            GSTATE.PRINT = false;
            if ~ischar(opt2)
                error('directory name must be a string');
            end
            if strcmp(opt2,'noprint')
                error('The output directory cannot be called ''noprint''');
            end
            GSTATE.DIR = opt2;
        else
            GSTATE.DIR = opt1;
            if ~strcmp(opt2,'noprint')
                error('Use ''noprint'' to avoid printing to file');
            end
            GSTATE.PRINT = false;
        end
end
  
stepf = 1/Nsymb;
GSTATE.FN = fftshift(-Nt/2:stepf:Nt/2-stepf);   % normalized frequencies
GSTATE.NSYMB = Nsymb;        % number of symbols
GSTATE.NT = Nt;              % points x symbol
GSTATE.NCH = Nch;            % number of channels
GSTATE.SYMBOLRATE = [];      % The symbolrate will be initialized in ELECTRICSOURCE

% The following fields will be initialized by CREATE_FIELD
GSTATE.FIELDX = [];          % x-component of the electric field. See CREATE_FIELD.
GSTATE.FIELDY = [];          % y-component of the electric field. See CREATE_FIELD.
GSTATE.FIELDX_TX =[];        % copy of GSTATE.FIELDX generated by CREATE_FIELD.
                             % Useful for back-to-back measurements.
GSTATE.FIELDY_TX = [];       % same as GSTATE.FIELDX_TX, but for y polarization.

GSTATE.DELAY = [];           % Overall delay cumulated in the optical line, normalized
                             % to 1/R. Size: [2,NCH].
GSTATE.DISP = [];            % cumulated dispersion [ps/nm] in the system.
                             % Size: [NPOL,NCH], being NPOL the number of 
                             % polarizations.
           
% The following fields will be initialized by LASERSOURCE
GSTATE.LAMBDA = [];          % Channels wavelength [nm]. Size: [1,NCH]
GSTATE.POWER = [];           % Transmitted Signal PEAK power [mW]. Size: [1,NCH]

if GSTATE.PRINT
    
    if ~exist(GSTATE.DIR,'dir');
       mkdir(GSTATE.DIR); % output directory       
    end
    if ~exist([GSTATE.DIR,'/',GSTATE.DIR,'.MOD/'],'dir');
       mkdir([GSTATE.DIR,'/',GSTATE.DIR,'.MOD/']);   % here print the power
    end
    if ~exist([GSTATE.DIR,'/',GSTATE.DIR,'.ANG/'],'dir');
       mkdir([GSTATE.DIR,'/',GSTATE.DIR,'.ANG/']);   % here print the phase
    end

    outfile = [GSTATE.DIR,'/simul_out'];   % matlab accepts both "/" and 
                                               % "\" under windows.

    nowt = clock;
    fid = fopen(outfile,'a');
    fprintf(fid,'++++++++++++++++++++++++++++++++++++++++\n');
    fprintf(fid,'++++       START OF SIMULATION      ++++\n');
    fprintf(fid,'++++                                ++++\n');
    if exist('OCTAVE_VERSION','var') == 1
        [pcname,t1] = system('hostname'); % under octave
    else
        [t1,pcname] = system('hostname'); % under matlab. Why? The answer 
    end                                   % is blowing in the wind...                  
    fprintf(fid,'++++ Hostname: %s',pcname);
    fprintf(fid,'++++ Date: %s %.2d:%.2d:%.2d\t\n',date,...
        nowt(4),nowt(5),round(nowt(6)));
    fprintf(fid,'++++++++++++++++++++++++++++++++++++++++\n\n\n');
    fprintf(fid,'========================================\n');
    fprintf(fid,'===             reset_all            ===\n');
    fprintf(fid,'========================================\n\n');
    fprintf(fid,'Global variable GSTATE initialized\n\n');
    fprintf(fid,'Nsymb = %6d\t (number of symbols)\n',Nsymb);
    fprintf(fid,'Nt   = %6d\t (points x symbol)\n',Nt);
    fprintf(fid,'Nch  = %6d\t (number of channels)\n',Nch);
    fprintf(fid,'Output directory = %s\n',GSTATE.DIR);
    fprintf(fid,'\n****************************************\n\n');

    fclose(fid);

    % Check for simul_out dimension

    simdir = dir(outfile);
    if (simdir.bytes > MAXBYTES)
        warning('optilux:reset_all',['The output file simul_out is very ',...
            'big: ',num2str(simdir.bytes),' bytes']);
        if nargout, wrn=true;end
    else
        if nargout, wrn=false;end
    end
end % end IF GSTATE.PRINT
