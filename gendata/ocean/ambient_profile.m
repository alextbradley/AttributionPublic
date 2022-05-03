function va = ambient_profile(pc, pw, z, vl, vu)
%Generate an two layer ambient profile with values vl and vu in the upper and lower layers repspectively. Returns profile at depth points z.  
va = zeros(length(pc), length(z));
for i = 1:length(pc) %loop over time points
    for j = 1:length(z) %loop over grid points

        va(i,j) = pointwise_ambient_profile(pc(i), pw(i), z(j), vl, vu);

    end
end
end

