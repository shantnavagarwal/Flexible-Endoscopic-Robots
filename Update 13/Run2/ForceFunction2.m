%% This function calculates force due to indentation; indentation velocity and friction.
function out = ForceFunction2(k, l, kv, u, x0, y0, ax, ay, avx, avy, j)
% k = stiffness; l = natural length of field (force applies if d > l); 
% kv is damping coefficient and u is coefficient of friction
% Beam B node = (x0, y0); force on node = (fx0, fy0)
% Beam A node   = (ax, ay); force on node = (afx, afy)
% beam A node velocity = (avx, avy)
d = sqrt((x0-ax)^2 + (y0-ay)^2);
% if d > l && (d < 2*l || ((j >= 1 && j <= 3) && ay <= y1 && ay <= y0 && ax < x0))
if d > l && d < 2*l
   % unit vector from node on A to point on line
    vectorx = (x0 - ax)/d;
    vectory = (y0 - ay)/d;
    % force magnitude
    f = -k*((d-l)^1);
    afx = -f*vectorx;
    afy = -f*vectory;
    fx0 = f*vectorx;
    fy0 = f*vectory;
    % velocity magnitude perpendicular to line
    vmag = avx * vectorx + avy * vectory;
    % damping force magnitude
    df = kv*vmag;
    afx = afx - df*vectorx;
    afy = afy - df*vectory;
    fx0 = fx0 + df*vectorx;
    fy0 = fy0 + df*vectory;
    region = 1;
%     disp([afx, afy])
%     out = [fx0, fy0, fx1, fy1, afx, afy]
else
    fx0 = 0;
    fy0 = 0;
    afx = 0;
    afy = 0;
    region = 0;
end
out = [fx0, fy0, afx, afy, region];
% disp(out)
end