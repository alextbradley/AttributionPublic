function va = pointwise_ambient_profile(pc, pw, z, vl, vu)
%return the ambient value at a depth z with the pycnocline centre at pc,
%width pw
d_low = pc - pw/2;
d_hi  = pc + pw/2;
if z< d_low
    va = vl;

elseif z > d_hi
    va = vu;

else
    va = vl + (vu - vl)/(d_hi - d_low) * (z - d_low);

end
end

