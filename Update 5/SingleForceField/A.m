%% Update 5
% Implements 1 circular forces
function [time,sig,f] = CIRC_POT(t,ne,le,e,ep,nx,lnp,x,xd)
f = [];
time = t;
nElem = 10;
n = 1 : 2 : (nElem*2+1);        % translational nodes, except the last (including last node causes strange behaviour)
ft(nElem+1,2) = 0;
f((nElem+1)*2, 3) = 0;
%% Circle 0
x0 = 8; %Circle center posn X
y0 = 0.75; % Circle center posn Y
r0 = 1; % Circle radius
k0 = .05; %Spring constant
for i = 1:nElem+1
    dist = (x(lnp(n(i),1))-x0)^2+(x(lnp(n(i),2))-y0)^2;
    if dist < r0^2
        ft(i,:) = ForceComponent(k0, x0, y0, x(lnp(n(i), 1)), x(lnp(n(i), 2)), r0);
    else
        ft(i, :) = [0, 0];
    end
    f(2*i-1, :) = [n(i) 1 ft(i,1)];
    f(2*i, :) = [n(i) 2 ft(i,2)];
%     disp("Circle 0")
%     disp(ft(i, :))
%     disp(f(2*i-1, :))
%     disp(f(2*i, :))
end
end