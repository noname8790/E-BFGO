function [Position,Fitness] = rand_devide(XB_Fitness,XB_position,K,n)
idx = randperm(K*n);
for k=1:K
    Position{k}=XB_position((k-1)*n+1:k*n,:); 
    Fitness{k}=XB_Fitness(1,(k-1)*n+1:k*n); 
end   
end

