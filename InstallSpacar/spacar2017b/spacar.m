function spacar(mode,filenm)
%SPACAR     Analysis of the dynamics of machines and mechanisms
%
%   SPACAR(MODE,FILENM) is the Matlab interface to the SPACAR program
%   that provides a tool for the analysis of the non-linear motion of
%   a mechanism as described in the input data file using a non-linear
%   finite element approach.
%
%   FILENM should specify the input file without the extension .DAT
%   Several output files will be created with the same name and different
%   extensions.
%
%   MODE selects one of the following
%        1 = Forward kinematic and dynamic analysis.
%        2 = Inverse dynamic analysis with setpoint generation.
%        3 = Linearization along predefined trajectory (of mode=2).
%        4 = Forward kinematic and dynamic analysis with linearization.
%        5 = Linearization for superposition (experimental).
%        6 = Linearization for superposition (experimental).
%   Using the negative values for MODE disables the plot of the 
%   mechanism during the simulation.
%
%   Consult the User Manual for more information.

% Note that this is only a help message as a MEX function will 
% actually be executed.

% Rev. 15-May-2001 R.Aarts.
