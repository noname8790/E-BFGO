function [Position,Fitness,lifetime] = devide(XB_Fitness,XB_position,XB_lifetime,K,n)
[Sorted_fitness,sorted_indexes]=sort(XB_Fitness);  
        for newindex=1:size(XB_position,1)
            Sorted_Position(newindex,:)=XB_position(sorted_indexes(newindex),:);
            Sorted_lifetime(newindex)=XB_lifetime(1,sorted_indexes(newindex));
        end
        for i=1:n
            XB_position(i,:)=Sorted_Position(6*i-5,:);
            XB_Fitness(1,i)=Sorted_fitness(1,6*i-5);
            XB_lifetime(1,i)=Sorted_lifetime(6*i-5);

            XB_position(i+n,:)=Sorted_Position(6*i-4,:);
            XB_Fitness(1,i+n)=Sorted_fitness(1,6*i-4);
            XB_lifetime(1,i+n)=Sorted_lifetime(6*i-4);

            XB_position(i+2*n,:)=Sorted_Position(6*i-3,:);
            XB_Fitness(1,i+2*n)=Sorted_fitness(1,6*i-3);
            XB_lifetime(1,i+2*n)=Sorted_lifetime(6*i-3);
            
            XB_position(i+3*n,:)=Sorted_Position(6*i-2,:);
            XB_Fitness(1,i+3*n)=Sorted_fitness(1,6*i-2);
            XB_lifetime(1,i+3*n)=Sorted_lifetime(6*i-2);

            XB_position(i+4*n,:)=Sorted_Position(6*i-1,:);
            XB_Fitness(1,i+4*n)=Sorted_fitness(1,6*i-1);
            XB_lifetime(1,i+4*n)=Sorted_lifetime(6*i-1);

            XB_position(i+5*n,:)=Sorted_Position(6*i,:);
            XB_Fitness(1,i+5*n)=Sorted_fitness(1,6*i);
            XB_lifetime(1,i+5*n)=Sorted_lifetime(6*i);
        end
        for k=1:K 
            Idx((k-1)*n+1:k*n)=k;
            Position{k}=XB_position(Idx==k,:); 
            Fitness{k}=XB_Fitness(1,Idx==k); 
            lifetime{k}=XB_lifetime(1,Idx==k);
        end  
end