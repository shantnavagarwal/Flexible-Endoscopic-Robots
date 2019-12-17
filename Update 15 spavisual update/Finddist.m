function out = Finddist(x0, y0, x1, y1, ax, ay)
% Beam B node 0 = (x0, y0)
% Beam B node 1 = (x1, y1)
% Beam A node = (ax, ay)
% Equation of line: ax + by + c = 0
eqa = y1 - y0;
eqb = x0 - x1;
eqc = -eqa*x0 - eqb*y0;
% Distance from center
d = abs((eqa*ax + eqb*ay + eqc)/sqrt(eqa^2 + eqb^2));
out = d;
end