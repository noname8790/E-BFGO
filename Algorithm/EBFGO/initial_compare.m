function [pp,pf,Gp,Gf,P,index,Elite] = initial_compare(Fitness,Position,K,Gp,Gf,pp,pf,n,index)
    for k=1:K
       [tf,I]= min(Fitness{k});
       if(tf<pf{k})
           pf{k}=tf;
           pp{k}=Position{k}(I,:);
           P(1,k)=0;
       else
           P(1,k)=1;
       end
    end
    for k=1:K
            if(pf{k}<Gf)
                r6=round(rand(1,1)*(K-1))+1;
                Fitness{r6}(1,n)=Gf;
                Position{r6}(n,:)=Gp;
                Gf=pf{k};
                Gp=pp{k};
            else
           Elite{index}=pp{k};
           index=index+1;
            end
    end  
end

