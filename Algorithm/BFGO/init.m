function pop = init(n_pop,d,Ub,Lb)
Z = zeros(n_pop, d);
% generate a random d-dimensional vector
Z(1, :) = rand(1, d);
% generate n_pop vectors using logistic 
for i=2:n_pop
    Z(i,:) = 4.0*Z(i-1,:).*(1-Z(i-1,:));
end
% map each component of z to the corresponding variable's value range
pop = zeros(n_pop, d);
for i=1:n_pop
    pop(i,:) = Lb + (Ub - Lb)*Z(i,:);
end
end

