function [Fx, Fy] = ForceComponent(k, x0, y0, x1, y1, r0)
d = sqrt((x1-x0)^2 + (y1-y0)^2);
dl = r0 - d;
Fx = ((x1-x0)/d)*k*dl;
Fy = ((y1-y0)/d)*k*dl;

