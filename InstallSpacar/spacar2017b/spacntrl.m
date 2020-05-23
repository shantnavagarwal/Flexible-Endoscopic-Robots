function spacntrl(filename,kpfun)
%SPACNTRL Includes the contribution of a controller in an LTV file.
%
%   SPACNTRL(FILENAME,KPFUN) read an LTV file, includes the contribution
%   of a controller and writes the results to a new LTV file.
%
%   FILENAME should specify the input file name WITHOUT the extension 
%   ".ltv". The output is written to a file FILENAME+"kp.ltv".
%
%   KPFUN is the name of an M-script that is called each time step to 
%   compute the control matrix gains. The calling syntax is 
%     KP = KPFUN(TIMESTEP)
%   where TIMESTEP is the number of the current time step and the size
%   of the return value KP has to match the number of actuators.
%   Before the first timestep the function is called with TIMESTEP=-1.
%
%   See also LTV, MRLTV.

% Note that this is only a help message as a MEX function will 
% actually be executed.

% Rev. 21-May-2001 R.Aarts.
