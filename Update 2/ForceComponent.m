function [Fx, Fy] = ForceComponent(k, x0, y0, x1, y1)
dx = x1 - x0;
dy = y1 - y0;
Fx = k*dx
Fy = k*dy
