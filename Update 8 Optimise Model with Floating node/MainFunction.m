clear
clc
%% Parameters to optimise
bendstiff = 0.0001; %in main_code.dat
longstiff = 300; %in main_code.dat
diststiff = 100; %in CIRC_POT.m
[B,L,D] = ndgrid(bendstiff, longstiff, diststiff);
%% Calculate Model force for Poke and Dist using force_analy
global poke
global dt
global ForceMat
global rodl
global BeamMat
rodl = 1;
dt = 0.5;
poke = .020;
dtsize = size(dt);
pokesize = size(poke);
ForceMat = zeros(dtsize(2), pokesize(2));
diam = 0.05;
thick = 0.005;
ymodul = 891640;
for d = 1:1:dtsize(2)
    dtin = dt(d);
    for p = 1:1:pokesize(2)
        ForceMat(d, p) = force_analy(ymodul, diam, thick, poke, dt);
    end
end
%% Grid search Matlab
fitresult = arrayfun(@(p1, p2, p3) Difference(p1, p2, p3), D, L, B);