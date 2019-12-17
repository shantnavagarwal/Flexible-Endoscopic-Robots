%% This function finds the node on which the doctor will apply force to insert the colonoscope into the large intestine.
% It does this by finding the node just below the entry point of the large
% testine. It also checks for the y value of each node of colonoscope to
% ensure correct node is being chosen
function node = findforcenod(lsa, lb)
% lsa: location of nodes of A from last point to enter LI to first point to
% enter LI
% nb : location of node of B should be 0, 0
mi = 100;
node = 0;
for i = 1:1:size(lsa, 1)
    n = lsa(i, :);
    if abs(n(1) - lb(1)) < 0.075 && lb(2) > n(2)
        if mi > lb(2) - n(2)
            node = i;
            mi = lb(2) - n(2);
        end
    end
end
node = 2*node - 1;
end