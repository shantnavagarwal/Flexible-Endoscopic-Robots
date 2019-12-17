function out = Finddist(x0, y0, x1, y1, ax, ay)
% Beam B node 0 = (x0, y0)
% Beam B node 1 = (x1, y1)
% Beam A node = (ax, ay)
% Equation of line: ax + by + c = 0
eqa = y1 - y0;
eqb = x0 - x1;
eqc = -eqa*x0 - eqb*y0;
linelength = sqrt((x1 - x0)^2 + (y1 - y0)^2);
% point on line where perendicular from node will intersect line
if eqa ~= 0 && eqb ~= 0
    bx = ((eqb^2)/eqa * ax - eqc - eqb*ay)/(eqa + (eqb^2)/ eqa);
    by = (-eqc - eqa*bx)/eqb;
else
    if eqa == 0
        bx = ax;
        by = -eqc/eqb;
    else
        bx = -eqc/eqa;
        by = ay;
    end  
end
bd0 = sqrt((bx - x0)^2 + (by - y0)^2);
bd1 = sqrt((bx - x1)^2 + (by - y1)^2);

% Distance from center
d = abs((eqa*ax + eqb*ay + eqc)/sqrt(eqa^2 + eqb^2));
% Check if in region
if bd0 <= linelength && bd1 <= linelength
    out = d;
else
    out = 200;
end
end