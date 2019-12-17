function repinsbf(filename,id,value,timestep)
%REPINSBF   Replace the value of a variable in a Spacar Binary data File
%
%   REPINSBF(FILENAME,ID,VALUE[,TIMESTEP]) replaces the value of one 
%   variable from the specified Spacar Binary data File (SBF).
%
%   FILENAME should specify the file name WITH the extension.
%   Note that the specified file is modified!!!!
%
%   ID is the name of the variable as e.g. is displayed with CHECKSBF.
%
%   VALUE is the new value of the variable. Its size should match
%   exactly to the original variable.
%
%   TIMESTEP is a required input parameter for time dependent data.
%   It is ignored for time independent data (i.e. header information).
%
%   See also CHECKSBF, GETFRSBF, LOADSBD, LOADSBM.

% Note that this is only a help message as a MEX function will 
% actually be executed.

% Rev. 21-May-2001 R.Aarts.
