function [good_position,poor_position] = fliped(Position,Fitness,K,Gp,Gf,pp,pf,n,Elite)
for k=1:K
    good_position{k}=Position{k}(1:n-2*n/5,:);
    poor_position{k}=Position{k}(n-2*n/5+1:n,:);
    F(1,:)=Gp;
    F(k+1,:)=pp{k};
end
%good
for k=1:K
    for i=1:n-n*2/5
        for d=1:size(Gp,2)
            r7=round(rand(1,1)*(K+1-1))+1;
            good_position{k}(i,d)=F(r7,d);
        end
    end
end
%poor
for k=1:K
    for i=1:n*2/5
        a=size(Elite,2);
     A=randperm(a,5);
     R= rand(1,5);
     R = R/sum(R);
     poor_position{k}(i,:)=R(1).*Elite{A(1)}+R(2).*Elite{A(2)}+R(3).*Elite{A(3)}+R(4).*Elite{A(4)}+R(5).*Elite{A(5)};
    end
end
end

