% Used to calculate coef of friction wrt sliding velocity
% coef of static friction is 1.2 times coef of sliding friction
% velocity threshold is 0.1
% Kinetic friction starts at 1.5 times vth
function uk = fricCoef(vlmag)
c1 = 1.2;
c2 = 1.5;
vlmag = abs(vlmag);
vth = 0.1;
if (vlmag < vth)
    uk = c1*vlmag/vth;
elseif (vlmag < c2*vth)
    vlmag = vlmag  - vth;
    uk = (c1-1)/(1 - c2)*(vlmag/vth)+1.2;
else
    uk = 1;
end    
end