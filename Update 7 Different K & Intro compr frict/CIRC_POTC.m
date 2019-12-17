%% Update 6
% Implements force field around a beam
function [time,sig,f] = CIRC_POT3(t,ne,le,e,ep,nx,lnp,x,xd)
    f = [];
    time = t;
    sig = [];

    nElem1 = 5;
    nElem2 = 5;
    n1 = 1 : 2 : (nElem1*2+1);        
    n2 = (nElem1*2+3):2:(nElem1*2+3)+nElem2*2;
    ft(nElem1+nElem2+2,2) = 0;
    f((nElem1+nElem2+2)*2, 3) = 0;
%     f((nElem1+1)*2, 3) = 0;
    
%%

    d = .95; %permissible distance
    k = .1;
    kv = 0.01;
    mu = 0.01;
    
    for i = 1:nElem1+1
        for j = 1:nElem2 
            f_out = ForceFunction(k, d, kv, mu, x(lnp(n2(j), 1)), x(lnp(n2(j), 2)), x(lnp(n2(j+1), 1)), x(lnp(n2(j+1), 2))...
                , x(lnp(n1(i), 1)), x(lnp(n1(i), 2)), xd(lnp(n1(i), 1)), xd(lnp(n1(i), 2)),3); 
            j_glob = nElem1+j+1;
            
            ft(i,:) = ft(i,:) + [f_out(5), f_out(6)];
            ft(j_glob,:) = ft(j_glob,:) + [f_out(1), f_out(2)];
            ft(j_glob+1,:) = ft(j_glob+1,:) + [f_out(3), f_out(4)];
            
%             f(2*j_glob-1, :) = [n2(j) 1 ft(j_glob,1)];
%             f(2*j_glob, :) = [n2(j) 2 ft(j_glob,2)];
            
        end
        f(2*i-1, :) = [n1(i) 1 ft(i,1)];
        f(2*i, :) = [n1(i) 2 ft(i,2)];
    end

    for j = 1:nElem2+1 %chk
        j_glob = nElem1+j+1;
        f(2*j_glob-1, :) = [n2(j) 1 ft(j_glob,1)];
        f(2*j_glob, :) = [n2(j) 2 ft(j_glob,2)];
    end
    %f
end
    