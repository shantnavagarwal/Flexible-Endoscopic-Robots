function loadsbd(filename)
%LOADSBD    Loads all data from a Spacar Binary Data (SBD) file.
%
%   LOADSBD(FILENAME) loads all data stored in the specified Spacar
%   Binary Data and creates or overwrites the accompanying workspace
%   variables: nx, ne, ndof, nddof, nkdof, nkddof, nxp, nep, 
%              lnp, le, ln, it, mode, rxyz, 
%              time, x, xd, xdd, e, ed, edd, fx, fxtot, fxgrav, sig, 
%              xcompl, xcompl, dx, de, dxc, dec, d2x, d2e
%
%   FILENAME should specify the file name WITHOUT the extension .SBD.
%
%   See also CHECKSBF, GETFRSBF, LOADSBM, REPINSBF.

% Note that this is only a help message as a MEX function will 
% actually be executed.

% Rev. 21-May-2001 R.Aarts.
