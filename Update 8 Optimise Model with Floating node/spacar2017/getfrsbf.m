function x=getfrsbf(filename,id,timesteps)
%GETFRSBF Extract a variable from a Spacar Binary data File
%
%   GETFRSBF(FILENAME,ID[,TIMESTEPS]) extracts one variable from
%   the specified Spacar Binary data File (SBF).
%
%   FILENAME should specify the file name WITH the extension.
%
%   ID is the name of the variable as e.g. is displayed with CHECKSBF.
%   To get the number of time steps take ID equal 'TDEF'.
%
%   TIMESTEPS is an optional input parameter for time dependent data.
%   It specifies the time steps at which the data should be read. If the
%   parameter is not specified, all time steps will be read.
%   If more than one time step is read, the output variable is a matrix 
%   of which the first index is the time step and all data are stored 
%   in a row.
%   This input parameter is ignored for time independent data (i.e. 
%   header information).
%
%   See also CHECKSBF, LOADSBD, LOADSBM, REPINSBF.

% Note that this is only a help message as a MEX function will 
% actually be executed.

% Rev. 16-Oct-2001 R.Aarts.
