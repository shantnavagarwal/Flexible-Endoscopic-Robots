global rodl
clear all
clc
global k
k = 1;
rodl = 1;
dt = 0.2:0.01:0.8;
poke = 0:.00001:.1;
dtsize = size(dt);
pokesize = size(poke);
ForceMat = zeros(dtsize(2), pokesize(2));
parfor d = 1:1:dtsize(2)
    dtin = dt(d);
    for p = 1:1:pokesize(2)
        ForceMat(d, p) = calcforce(dtin, poke(p));
    end
end
        