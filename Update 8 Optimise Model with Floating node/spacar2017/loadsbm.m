function loadsbm(filename)
%LOADSBM    Loads all data from a Spacar Binary Matrix data (SBM) file.
%
%   LOADSBM(FILENAME) loads all data stored in the specified Spacar
%   Binary Matrix data and creates or overwrites the accompanying 
%   workspace variables: nddof, ndof, nkddof, nkdof, nnom, 
%                        time, m0, c0, d0, k0, n0, g0, ak0, bk0, b0.
%
%   FILENAME should specify the file name WITHOUT the extension .SBM.
%
%   Note that in the current SPACAR release the structural stiffness
%   matrix K0 and the dynamic stiffening matrix N0 are added in k0.
%
%   See also CHECKSBF, GETFRSBF, LOADSBD, REPINSBF.

% Note that this is only a help message as a MEX function will 
% actually be executed.

% Rev. 21-May-2001 R.Aarts.
