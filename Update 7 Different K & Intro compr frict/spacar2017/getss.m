function [A,B,C,D] = getss(filename, ts)
%GETSS    Loads system matrices from a Spacar LTV file.
%
%   [A,B,C,D] = GETSS(FILENAME, TS) loads the matrices from the specified
%   Spacar LTV file at the specified time step. The output is stored in 
%   four variables. 
%
%   SYS = GETSS(FILENAME, TS) loads the matrices from the specified
%   Spacar LTV file at the specified time step. The output is a single 
%   SS object generated from these data. 
%
%   FILENAME should specify the file name WITHOUT the extension .LTV.
%
%   When no time step TS is specified the first is assumed. 
%
%   See also CHECKSBF.

% Rev. 26-Jan-2007 R.Aarts.

error(nargchk(1,2,nargin));
error(nargchk(0,4,nargout));
if nargout==2 || nargout==3
    error('Output is either one SS object or 4 state space matrices')
end

if nargin<2
    ts = 1;
else
    ts = ts(1,1);
end

if ~ischar(filename)
   error('First parameter must be filename.')
end

if length(filename)<=4 || ~strcmpi(filename(end-3:end), '.ltv')
    filename = [ filename '.ltv' ];
end

if ~exist(filename, 'file')
    error(['SPACAR file ' filename ' not found!']);
end
        
A = getfrsbf(filename,'A',ts);
B = getfrsbf(filename,'B',ts);
C = getfrsbf(filename,'C',ts);
if getfrsbf(filename,'DFT')
    D = getfrsbf(filename,'D',ts);
else
    D = [];
end

if nargout<2
    A=ss(A,B,C,D);
    clear B C D
end

