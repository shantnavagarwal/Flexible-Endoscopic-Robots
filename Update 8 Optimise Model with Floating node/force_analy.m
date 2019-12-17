function f_y=force_analy(E,D,thick,y_pen,dist)
% E: Young's Modulus, D: outer Diameter
% Y_pen: peneteration, dist: Distance from one end of line segment
d=D-thick*2;
I=pi*(D^4-d^4)/64;
f_y=abs(y_pen*6*E*I/(dist*(1-dist)*(2*dist^2-2*dist)));
end