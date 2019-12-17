f = [];
time = 5;
nElem = 20;
n = 1 : 2 : (nElem*2-1);        % translational nodes, except the last
ft(nElem,2) = 0;
x0 = 1; %Circle center posn X
y0 = 1; % Circle center posn Y
r = 2; % Circle radius
k = 1; %Spring constant

for i = 1:nElem
    dist = x(lnp(n(i),1)-x0)^2+(x(lnp(n(i),2)-y0))^2;
    if dist < r^2
        ft(i,:) = sqrt(dist)*k;
    else
        ft(i,:) = 0;
    end
end
       