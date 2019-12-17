function [time,sig,f] = CIRC_POT(t,ne,le,e,ep,nx,lnp,x,xd)

f = [];
time = t;
nElem = 10;
%inp_for = 0;
n = 1 : 2 : (nElem*2-1);        % translational nodes, except the last
ft(nElem,2) = 0;
x0 = 3; %Circle center posn X
y0 = 1; % Circle center posn Y
r = 2; %Circle radius
k = .01; %Spring constant

% if (x(lnp(1,1))-x0)^2+(x(lnp(1,2))-y0)^2 >= r^2
%     f = [f; [1 1 inp_for]];
%     %f = [f; [1 2 0]];
% else
%     ft(1,: ) = ForceComponent(k, x0, y0, x(lnp(1, 1)), x(lnp(1, 2)));
%     f = [f; [1 1 inp_for+ft(1,1)]];
%     %f = [f; [1 2 ft(1,2)]];
% end
    
for i = 1:nElem
    %lnp(n(i),1)
    dist = (x(lnp(n(i),1))-x0)^2+(x(lnp(n(i),2))-y0)^2;
    if dist < r^2
        ft(i,:) = ForceComponent(k, x0, y0, x(lnp(n(i), 1)), x(lnp(n(i), 2)));
    else
        ft(i, :) = [0, 0];
    end
%     if i == 1
%         ft(i, 1) = ft(i, 1)+inp_for;
%     end
    f = [f; [n(i) 1 ft(i,1)]];
    f = [f; [n(i) 2 ft(i,2)]];
end
% if flag
%     f(1,3)= f(1,3)+inp_for;

sig = [];
end
        
        
        
        