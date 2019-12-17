function BezierFit(Q, n)
% BezierFitDemo.m  Fitting Data using Picecewire G1 Cubic Bézier Curves
%
% Demonstration of MATLAB code from:
% Lane, Edward J. Fitting Data Using Piecewise G1 Cubic Bezier Curves.
% Thesis, NAVAL POSTGRADUATE SCHOOL MONTEREY CA, 1995.
%
% This code depends on the geom2D Toolbox which requires MATLAB R2014b or
% later.
% https://www.mathworks.com/matlabcentral/fileexchange/7844-geom2d


% Demonstrate two different Bezier curve fits
% Array Q is the set of data points to which we want to fit
% a piecewise G1 continuous cubic Bezier curve
% (Replace Q with your data.  You will need to adjust n = number of knots)

% Now try to do a piecewise cubic Bézier fit to Q starting with n knot points
% This computes the globally optimized only (GOO) curvve.
% (The plotting of the IG curve in iguess0.m can be commented out.)
Qt = Q';
[IG, k] = iguess0(Qt, n);

% Improve the fit: Segmentally Optimum Only Curve (SOO)
SOC = segop(k, Qt, IG);

% Improve it again: Segmentally then Globally Optimized Curve (SGO)
GOC = globop(SOC, Qt, 0, k);

% plot SGO curve using internal routine
figure;
hp = poplt(GOC, Qt);
title('Plot of SGO curve')

% This section of code simply repeats the plot above (poplt) in a
% manner that makes it easier to see the piecewise cubic Beizier sections

% Get the Bézier control points of the curve fit
Cnew = cpoints(GOC)';
P    = knots(Qt, k)';

% Plot fitted, segment by segment using Geom2D toolkit
figure; 
hold on;
for i = 1:3:length(Cnew)-3
    disp(i:i+3)
    drawBezierCurve(Cnew(i:i+3,:));      % Fitted cubic bezier segment
end
hc = plot(Cnew(:,1), Cnew(:,2), 'o-');   % Control points
ho = plot(Q(:,1), Q(:,2), 'k.');           % Original data
hn = plot(P(:,1), P(:,2), 'kx', 'markerSize', 10);   % Original knot points
legend([ho hn hc], 'Original data', ...
    sprintf('Original guess n = %d knots',n), ...
    'New control points');
set(gcf, 'color', 'w');