clear
clc
tic
parfor i = 1:4
    if (mod(i, 4) == 1)
        disp('1')
        disp('Constant')
        spacar(1, 'main_code1')
    elseif (mod(i, 4) == 2)
        disp('2')
        disp('Exponential')
        spacar(1, 'main_code2')
    elseif (mod(i, 4) == 3)
        disp('3')
        disp('Square')
        spacar(1, 'main_code3')
    elseif (mod(i, 4) == 4)
        disp('4')
        disp('logarithmic')
        spacar(1, 'main_code4')
    end
end
toc