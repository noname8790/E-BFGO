function [St]=Checkbound(St,Stlow_bound,Stup_bound,Np,Dimension,G)

    for i=1:Dimension
        for k=1:Np
            if St(k,i)<Stlow_bound
                St(k,i)=rand*(Stup_bound-Stlow_bound)+Stlow_bound;%*(1+G/10);
            else if St(k,i)>Stup_bound
                St(k,i)=Stup_bound-rand*(Stup_bound-Stlow_bound);%/(1+G/10);
                end
            end
        end
    end
end