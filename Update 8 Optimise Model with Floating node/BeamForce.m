function F_y = BeamForce(k,lon_stiff,ben_stiff)
%     global tip_forcex
    global tip_forcey
    global k_val
    global poke
    global dt
    k_val=k;
    %foo=fopen('main_code.dat','w');
    %save main_code.dat x -ascii
    dtsize = size(dt);
    pokesize = size(poke);
%     F_x = zeros(dtsize(2), pokesize(2));
    F_y = zeros(dtsize(2), pokesize(2));
    for i =1:1:dtsize(2)
        dist = dt(i);
%         f_x = zeros(1, pokesize);
        f_y = zeros(pokesize);
        for j = 1:1:pokesize(2)
            depth = poke(j);    
            dat_writer(depth,dist,lon_stiff,ben_stiff);
            spacar(1,'main_code');
%             f_x(1, j)=tip_forcex;
            f_y(1, j)=tip_forcey;
        end
%         F_x(i, :) = f_x;
        F_y(i, :) = f_y;
    end
end