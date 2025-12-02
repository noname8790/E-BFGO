function mag = norm1(v)
sv = v.* v;     
dp = sum(sv);    
mag = sqrt(dp); 
end

