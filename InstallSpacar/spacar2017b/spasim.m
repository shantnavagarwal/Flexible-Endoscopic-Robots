function [sys,x0,str,ts]=spasim(t,x,u,flag,filename,plot)
%SPASIM     Simulation of the dynamics of machines and mechanisms
%
%   [SYS,X0,STR,TS]=SPASIM(T,X,U,FLAG,FILENAME,PLOT) is the Simulink 
%   interface to the SPACAR program that provides a tool for the
%   simulation of the non-linear motion of a mechanism as described
%   in the input data file using a non-linear finite element approach.
%
%   Consult the User Manual for more information.
%
%   Most parameters have their common meaning for S-functions.
%        
%   FILENM is the name of the file to read (without the extension DAT).
%   PLOT   is switch to disable/enable plotting of the mechanism

% Note that this is only a help message as a MEX function will 
% actually be executed.

% Rev. 27-Sep-2002 R.Aarts.
