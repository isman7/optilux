function ber=q2ber(Q)

%BER=Q2BER(Q) Convert the Q-factor in bit-error rate
%	BER=Q2BER(Q) Converts the Q-factor [dB] in bit-error rate 
%	using the formula:
%
%   BER = 0.5 * erfc( Qlin / sqrt(2) )
%
%   where Qlin=10^(Q/20).
%
%	As a reference (BER -> Q [dB}): 
%						1e-3 -> 9.80		8e-4 ~= 10
%						1e-5 -> 12.60		2e-4 ~= 11
%						1e-9 -> 15.56		
%
%   Note: ber2q(q2ber(Q)) = Q
%
%   See also BER2Q
%
%   Author: Paolo Serena, 2010
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

ber=0.5*erfc( 10.^(Q/20) / sqrt(2));