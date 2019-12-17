%% This function calculates force due to indentation; indentation velocity and friction.
function out = ForceFunction(k, l, kv, u, x0, y0, x1, y1, ax, ay, avx, avy)
% k = stiffness; l = natural length of field (force applies if d > l); 
% kv is damping coefficient and u is coefficient of friction
% Beam B node 0 = (x0, y0); force on node = (fx0, fy0)
% Beam B node 1 = (x1, y1); force on node = (fx1, fy1)
% Beam A node   = (ax, ay); force on node = (afx, afy)
% beam A node velocity = (avx, avy)
% Equation of line: ax + by + c = 0
eqa = y1 - y0;
eqb = x0 - x1;
eqc = -eqa*x0 - eqb*y0;
inregion = 0;
% length of line segment
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
if bd0 <= linelength && bd1 <= linelength && d < 2*l
    inregion = 1;
end
if bd0 <= linelength && bd1 <= linelength && d > l && d < 2*l
%     disp('gussa yay')
   % unit vector from node on A to point on line
    vectorx = (bx - ax)/d;
    vectory = (by - ay)/d;
    % force magnitude
    f = -k*((d-l)^1.5);
    afx = -f*vectorx;
    afy = -f*vectory;
    fx0 = bd1/linelength*f*vectorx;
    fy0 = bd1/linelength*f*vectory;
    fx1 = bd0/linelength*f*vectorx;
    fy1 = bd0/linelength*f*vectory;
    % velocity magnitude perpendicular to line
    vmag = avx * vectorx + avy * vectory;
    % damping force magnitude
    df = kv*vmag;
    afx = afx - df*vectorx;
    afy = afy - df*vectory;
    fx0 = fx0 + bd1/linelength*df*vectorx;
    fy0 = fy0 + bd1/linelength*df*vectory;
    fx1 = fx1 + bd0/linelength*df*vectorx;
    fy1 = fy1 + bd0/linelength*df*vectory;
    % friction
    % unit vector along line
    linex = (x1 - x0)/linelength;
    liney = (y1 - y0)/linelength;
    % velocity magnitude along line
    vlmag = avx*linex + avy*liney;
    % friction force
    uk = fricCoef(vlmag);
    fricf = -sign(vlmag)*uk*u*sqrt(afx^2+afy^2);
    afx = afx + fricf*linex;
    afy = afy + fricf*liney;
    out = [fx0, fy0, fx1, fy1, afx, afy, inregion];
%     disp(out)
else
    fx0 = 0;
    fx1 = 0;
    fy0 = 0;
    fy1 = 0;
    afx = 0;
    afy = 0;
end
out = [fx0, fy0, fx1, fy1, afx, afy, inregion];
% disp(out)
end