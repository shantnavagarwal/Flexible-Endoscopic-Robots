function [sys,x0,str,ts]=ltv(t,x,u,flag,filename[,mode])
%LTV        Read data from an LTV file for use in a Simulink model
%
%   [SYS,X0,STR,TS]=LTV(T,X,U,FLAG,FILENAME[,MODE]) is the Simulink 
%   interface to read the Linear Time-Varying data from an LTV file
%   for use in a Simulink model. Its behaviour depends strongly on the
%   (optional) parameter MODE.
%
%   Consult the User Manual for more information.
%
%   Most parameters have their common meaning for S-functions.
%        
%   FILENAME is the name of the file to read (without the extension LTV).
%
%   MODE is an optional one dimensional array:
%     IMODE = MODE(1): Mode (behaviour) of the Simulink block:
%        0: Represent an LTV block, that means use the state space 
%           matrices in the LTV file to simulate a Linear Time-Varying
%           system. Depending on the input file the system is an 
%           "ordinary" state space system, or it is the result from
%           modal analysis.
%           In case the parameter MODE in not specified, this is the
%           default behaviour.
%       11: Read vector U0 from the LTV file. 
%       12: Read vector Y0 from the LTV file. 
%       13: Read vector SIG0 from the LTV file. 
%       21: Read matrix M0 from the LTV file. The output in Simulink
%           is the matrix in a linear representation.
%       31: Multiply the input vector with matrix M0. In case the size
%           of the input vector is smaller than M0, MODE(3) should
%           specify the size.
%     INTPOL = MODE(2): Interpolation type of the output:
%        0: Stepwise, i.e. constant around the time steps found in the 
%           LTV file.
%        1: Linear (this is the default).
%        2: Cubic spline (not available for IMODE==0).
%     NMODES = MODE(3): Number of modes. This parameters is meaningful
%       for IMODE==0 and the LTV file is suitable for modal reduction,
%       or for IMODE==31 (see above).
%       <0: All modes.
%      >=0: (Reduced) number of modes to be used.
%
%   See also MRLTV, SPACNTRL.

% Note that this is only a help message as a MEX function will 
% actually be executed.

% Rev. 21-May-2001 R.Aarts.
