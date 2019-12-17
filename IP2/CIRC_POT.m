%% Update 14
% Implements force field around a beam
function [time,sig,f] = CIRC_POT(t,ne,le,e,ep,nx,lnp,x,xd)
    global ff %debug
    f = [];
    ff = [];
    time = t;
    sig = [];
    %% Parameters
    d = .025; %permissible distance
    k = 150;
    kv = 0.01;
    mu = 0.2;
    nElem1 = 10;
    nElem2 = 23;
    %% Not to be changed
    n1 = 1 : 2 : (nElem1*2+1);        
    n2 = (nElem1*2+3):2:(nElem1*2+3)+nElem2*2;
    ft(nElem1+nElem2+2,2) = 0;
    f((nElem1+nElem2+2)*2, 3) = 0;
%     f((nElem1+1)*2, 3) = 0;
    regionA = zeros(nElem1+1,1);
    
    
    for i = 1:nElem1+1
        min_elem=100; %min_dist from any element
        min_pos_e = -1;
%         for_from_nearest=0;
%         bool_force = 0; %whether force is applied at all
        for j = 1:nElem2
            
            dist = Finddist(x(lnp(n2(j), 1)), x(lnp(n2(j), 2)), x(lnp(n2(j+1), 1)), x(lnp(n2(j+1), 2))...
                , x(lnp(n1(i), 1)), x(lnp(n1(i), 2))); %finds distance from each element
            
            f_out = ForceFunction(k, d, kv, mu, x(lnp(n2(j), 1)), x(lnp(n2(j), 2)), x(lnp(n2(j+1), 1)), x(lnp(n2(j+1), 2))...
                , x(lnp(n1(i), 1)), x(lnp(n1(i), 2)), xd(lnp(n1(i), 1)), xd(lnp(n1(i), 2)));
            if dist > d && dist < min_elem
                min_elem = dist;
                min_pos_e = j;
            end
%             if f_out ~= 0
%                 disp([i j])
%                 disp(t)
%                 disp(f_out)
%             end
%             if dist < 2*d
            j_glob = nElem1+j+1;
%             disp(f_out)
            ft(i,:) = ft(i,:) + [f_out(5), f_out(6)];
            ft(j_glob,:) = ft(j_glob,:) + [f_out(1), f_out(2)];
            ft(j_glob+1,:) = ft(j_glob+1,:) + [f_out(3), f_out(4)];
            if f_out(7)==1
                regionA(i)=1;
            end
%         end
%             f(2*j_glob-1, :) = [n2(j) 1 ft(j_glob,1)];
%             f(2*j_glob, :) = [n2(j) 2 ft(j_glob,2)];            
        end
        
%         f(2*i-1, :) = [n1(i) 1 ft(i,1)];
%         f(2*i, :) = [n1(i) 2 ft(i,2)];
        
        if regionA(i)==0 && (x(lnp(n1(i), 1))>=0.1 || x(lnp(n1(i), 2))>=0.0)
            min_node=Inf; %min_dist from any node
            min_pos_n=0;
            for j=1:nElem2+1
                cur_dis=calc_dist(x(lnp(n2(j), 1)), x(lnp(n2(j), 2)), x(lnp(n1(i), 1)), x(lnp(n1(i), 2)));
                if cur_dis<min_node
                    min_node=cur_dis;
                    min_pos_n=j;
                end
            end
            if min_node>d && min_node<2*d %if min_dist from the node
                f_out2=ForceFunction2(k, d, kv, mu, x(lnp(n2(min_pos_n), 1)), x(lnp(n2(min_pos_n), 2))...
                    , x(lnp(n1(i), 1)), x(lnp(n1(i), 2)), xd(lnp(n1(i), 1)), xd(lnp(n1(i), 2)));
                j_glob = nElem1+min_pos_n+1;
                ft(i,:) = ft(i,:) + [f_out2(3), f_out2(4)];
                ft(j_glob,:) = ft(j_glob,:) + [f_out2(1), f_out2(2)];
            
            elseif min_pos_e>-1 && min_elem < min_node %if no dist within given limit, then force preferably from element
                j=min_pos_e;
                f_out = ForceFunction3(k, d, kv, mu, x(lnp(n2(j), 1)), x(lnp(n2(j), 2)), x(lnp(n2(j+1), 1)), x(lnp(n2(j+1), 2))...
                , x(lnp(n1(i), 1)), x(lnp(n1(i), 2)), xd(lnp(n1(i), 1)), xd(lnp(n1(i), 2)));
%                 f_out = for_from_nearest;
                j_glob = nElem1+j+1;
                ft(i,:) = ft(i,:) + [f_out(5), f_out(6)];
                ft(j_glob,:) = ft(j_glob,:) + [f_out(1), f_out(2)];
                ft(j_glob+1,:) = ft(j_glob+1,:) + [f_out(3), f_out(4)];
            
            elseif min_node > d %if not in any of the regions, force applicable from the nearest node
                f_out2=ForceFunction2(k, d, kv, mu, x(lnp(n2(min_pos_n), 1)), x(lnp(n2(min_pos_n), 2))...
                    , x(lnp(n1(i), 1)), x(lnp(n1(i), 2)), xd(lnp(n1(i), 1)), xd(lnp(n1(i), 2)));
                j_glob = nElem1+min_pos_n+1;
                ft(i,:) = ft(i,:) + [f_out2(3), f_out2(4)];
                ft(j_glob,:) = ft(j_glob,:) + [f_out2(1), f_out2(2)];
            end
        
        end
        f(2*i-1, :) = [n1(i) 1 ft(i,1)];
        f(2*i, :) = [n1(i) 2 ft(i,2)];
    end

    for j = 1:nElem2+1 %chk
        j_glob = nElem1+j+1;
        f(2*j_glob-1, :) = [n2(j) 1 ft(j_glob,1)];
        f(2*j_glob, :) = [n2(j) 2 ft(j_glob,2)];
%         if f(2*j_glob,3)>1 || f(2*j_glob-1,3)>1
%             f(2*j_glob,3)
%         end
    end
    %f
%     disp(t)
%     disp(f_out)
%     disp(f)
%     ff = [ff; f];
end
    