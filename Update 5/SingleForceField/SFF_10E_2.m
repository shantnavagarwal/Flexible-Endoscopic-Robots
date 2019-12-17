%% Update 4
% Implements 2 circular forces
function [time,sig,f] = CIRC_POT(t,ne,le,e,ep,nx,lnp,x,xd)
f = [];
time = t;
nElem = 100;
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
%% Circle 1
x1 = 7.5; %Circle center posn X
y1 = -.7; % Circle center posn Y
r1 = 0.45; % Circle radius
k1 = .05; %Spring constant
for i = 1:nElem+1
    dist = (x(lnp(n(i),1))-x1)^2+(x(lnp(n(i),2))-y1)^2;
    if dist < r1^2
        ft(i,:) = ForceComponent(k1, x1, y1, x(lnp(n(i), 1)), x(lnp(n(i), 2)), r1);
    else
        ft(i, :) = [0, 0];
    end
    f(2*i-1, 3) = f(2*i-1, 3) + ft(i, 1);
    f(2*i, 3) = f(2*i, 3) + ft(i, 2);
%     disp("Circle 1")
%     disp(ft(i, :))
%     disp(f(2*i-1, :))
%     disp(f(2*i, :))
end
%% Circle 2
x2 = 8.5; %Circle center posn X
y2 = -.95; % Circle center posn Y
r2 = 0.5; % Circle radius
k2 = .08; %Spring constant
for i = 1:nElem+1
    dist = (x(lnp(n(i),1))-x2)^2+(x(lnp(n(i),2))-y2)^2;
    if dist < r2^2
        ft(i,:) = ForceComponent(k2, x2, y2, x(lnp(n(i), 1)), x(lnp(n(i), 2)), r2);
    else
        ft(i, :) = [0, 0];
    end
    f(2*i-1, 3) = f(2*i-1, 3) + ft(i, 1);
    f(2*i, 3) = f(2*i, 3) + ft(i, 2);
%     disp("Circle 2")
%     disp(ft(i, :))
%     disp(f(2*i-1, :))
%     disp(f(2*i, :))
end
sig = [];
end
        
        
        
        