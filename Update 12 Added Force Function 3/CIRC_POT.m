%% Update 6
% Implements force field around a beam
function [time,sig,f] = CIRC_POT(t,ne,le,e,ep,nx,lnp,x,xd)
    global ff %debug
    f = [];
    ff = [];
    time = t;
    sig = [];
    %% Parameters
    d = .025; %permissible distance
    k = 200;
    kv = 0.01;
    mu = 0.01;
    nElem1 = 20;
    nElem2 = 23;
    %% Not to be changed
    n1 = 1 : 2 : (nElem1*2+1);        
    n2 = (nElem1*2+3):2:(nElem1*2+3)+nElem2*2;
    ft(nElem1+nElem2+2,2) = 0;
    f((nElem1+nElem2+2)*2, 3) = 0;
%     f((nElem1+1)*2, 3) = 0;
    regionA = zeros(nElem1+1,1);
    
    
    for i = 1:nElem1+1
        for j = 1:nElem2
%             if i == 6 && j == 2 && t >1
%                 disp('fdas')
%             end
          
            f_out = ForceFunction(k, d, kv, mu, x(lnp(n2(j), 1)), x(lnp(n2(j), 2)), x(lnp(n2(j+1), 1)), x(lnp(n2(j+1), 2))...
                , x(lnp(n1(i), 1)), x(lnp(n1(i), 2)), xd(lnp(n1(i), 1)), xd(lnp(n1(i), 2)));
%             if f_out ~= 0
%                 disp([i j])
%                 disp(t)
%                 disp(f_out)
%             end
            j_glob = nElem1+j+1;
%             disp(f_out)
            ft(i,:) = ft(i,:) + [f_out(5), f_out(6)];
            ft(j_glob,:) = ft(j_glob,:) + [f_out(1), f_out(2)];
            ft(j_glob+1,:) = ft(j_glob+1,:) + [f_out(3), f_out(4)];
            if f_out(7)==1
                regionA(i)=1;
            end
%             f(2*j_glob-1, :) = [n2(j) 1 ft(j_glob,1)];
%             f(2*j_glob, :) = [n2(j) 2 ft(j_glob,2)];
            
        end
%         f(2*i-1, :) = [n1(i) 1 ft(i,1)];
%         f(2*i, :) = [n1(i) 2 ft(i,2)];
        if regionA(i)==0
            min=Inf;
            min_pos=0;
            for j=1:nElem2+1
                cur_dis=calc_dist(x(lnp(n2(j), 1)), x(lnp(n2(j), 2)), x(lnp(n1(i), 1)), x(lnp(n1(i), 2)));
                if cur_dis<min
                    min=cur_dis;
                    min_pos=j;
                end
            end
            f_out2=ForceFunction2(k, d, kv, mu, x(lnp(n2(min_pos), 1)), x(lnp(n2(min_pos), 2))...
                , x(lnp(n1(i), 1)), x(lnp(n1(i), 2)), xd(lnp(n1(i), 1)), xd(lnp(n1(i), 2)));
            j_glob = nElem1+min_pos+1;
            ft(i,:) = ft(i,:) + [f_out2(3), f_out2(4)];
            ft(j_glob,:) = ft(j_glob,:) + [f_out2(1), f_out2(2)];
        end
        %force to keep endoscope straight
        if x(lnp(n1(i), 2))<0 && x(lnp(n1(i), 1))<=-0.01
           ft(i,:) = ft(i,:) + [ForceFunction3(x(lnp(n1(i), 1))), 0]; 
        end
        if x(lnp(n1(i), 2))<0 && x(lnp(n1(i), 1))>=0.01
           ft(i,:) = ft(i,:) + [ForceFunction3(x(lnp(n1(i), 1))), 0]; 
        end
        f(2*i-1, :) = [n1(i) 1 ft(i,1)];
        f(2*i, :) = [n1(i) 2 ft(i,2)];
    end

    for j = 1:nElem2+1 %chk
        j_glob = nElem1+j+1;
        f(2*j_glob-1, :) = [n2(j) 1 ft(j_glob,1)];
        f(2*j_glob, :) = [n2(j) 2 ft(j_glob,2)];
        if f(2*j_glob,3)>1 || f(2*j_glob-1,3)>1
            f(2*j_glob,3)
        end
    end
    %f
%     disp(t)
%     disp(f_out)
%     disp(f)
%     ff = [ff; f];
end
    