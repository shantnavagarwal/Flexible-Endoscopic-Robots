%% Deprecated
function force = calcforce(d, p)
disp('calcforce.m is Deprecated. See force_analy.m')
global rodl
rodl = 1;
k = 1;
l1 = sqrt(d^2 + p^2);
l2 = sqrt((rodl-d)^2 + p^2);
force = k*(l1-d)*(p/l1) + k*(l2 - (rodl-d))*(p/l2);
end