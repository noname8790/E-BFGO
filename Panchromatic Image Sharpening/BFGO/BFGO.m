function [Gf,Gp,Convergence_curve] = BFGO(func,D,N,T,lb,ub)
%% initialize bamboo population
rand('state',sum(100*clock)*rand(1));
XB_position=init(N,D,ub,lb);
for i=1:N
    XB_fitness(1,i)=func(XB_position(i,:));
end
n=5;
K=N/n;% number of clusters
Gf=Inf;
Gp=zeros(1,D);
pp=cell(1,K);
pf=cell(1,K);
for k=1:K
    pf{k}=Inf;
end
index=1;
[Position,Fitness]=devide(XB_fitness,XB_position,K,n);% divided into k groups
[pp,pf,Gp,Gf,P,index,Elite]=initial_compare(Fitness,Position,K,Gp,Gf,pp,pf,n,index);% update the best individual
Convergence_curve(1)=Gf;
t=2;
%% iteration
while t<T+1
     if(t>10&&Convergence_curve(t-1)==Convergence_curve(t-3))
        [good_position,poor_position]=fliped1(Position,Fitness,K,Gp,Gf,pp,pf,n,Elite);
        for k=1:K
            for i=1:3*n/5
                Flag4ub=good_position{k}(i,:)>ub;
                Flag4lb=good_position{k}(i,:)<lb;
                good_position{k}(i,:)=(good_position{k}(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;  
                good_fitness{k}(1,i)=func(good_position{k}(i,:));
            end
            for i=1:2*n/5
                Flag4ub=poor_position{k}(i,:)>ub;
                Flag4lb=poor_position{k}(i,:)<lb;
                poor_position{k}(i,:)=(poor_position{k}(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
                poor_fitness{k}(1,i)=func(poor_position{k}(i,:));
            end
            for i=1:3*n/5
            if(good_fitness{k}(1,i)<Fitness{k}(1,i))
                Position{k}(i,:)=good_position{k}(i,:);
                Fitness{k}(1,i)=good_fitness{k}(1,i);
            end
            end
            Position{k}(n-2*n/5+1:n,:)=poor_position{k};
            Fitness{k}(1,n-2*n/5+1:n)=poor_fitness{k};
        end
         [pp,pf,Gp,Gf,P,index,Elite]=compare(Fitness,Position,K,Gp,Gf,pp,pf,n,index,Elite);% update the best individual
     end   
   
    
    
    if(sum(P)==K)
        [Position,Fitness]=rand_devide(XB_fitness,XB_position,K,n);% randomly divided into k groups
    end
%% bamboo whip extension stage
    a=1-t*((1)/T);b=1+4*t/T;
        for k=1:K
            C{k}=zeros(1,D);
            for i=1:n
            C{k}=C{k}+Position{k}(i,:);
            end
            C{k}=C{k}./n;
        for i=1:n
            r1=rand();
            if(r1<0.4)
               cossita=dot(Gp,Position{k}(i,:))/norm1(Gp)/norm1(Position{k}(i,:));
               Position{k}(i,:)=Gp+2*a*(rand()*Gp-Position{k}(i,:))*cossita;
            else if(r1>0.4&&r1<0.8) 
                    cossita=dot(pp{k},Position{k}(i,:))/norm1(pp{k})/norm1(Position{k}(i,:));
                    Position{k}(i,:)=Gp+2*a*(rand()*pp{k}-Position{k}(i,:))*cossita;
                else
                    cossita=dot(C{k},Position{k}(i,:))/norm1(C{k})/norm1(Position{k}(i,:));
                    Position{k}(i,:)=C{k}+2*a*(rand()*C{k}-Position{k}(i,:))*cossita;
                end
            end
            Flag4ub=Position{k}(i,:)>ub;
            Flag4lb=Position{k}(i,:)<lb;
            Position{k}(i,:)=(Position{k}(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;   
            Fitness{k}(1,i)=func(Position{k}(i,:));
        end  
        end
        [pp,pf,Gp,Gf,P,index,Elite]=compare(Fitness,Position,K,Gp,Gf,pp,pf,n,index,Elite);% update the best individual 
    
%% bamboo shoot growth stage
    b=2*rand();
    sita=2;
    for k=1:K
        XB(k).x(t,:)=Gp*exp(b/sita*t^sita);
        XB(k).dertx=XB(k).x(t,:)-XB(k).x(t-1,:);
    end   
    for k=1:K
        for i=1:n
            r3=round(rand(1,1)*(n-1))+1;
            XB(k).cx(i,:)=XB(k).dertx./(Gp-Position{k}(i,:));
            XB(k).dis(i,:)=1-abs((Position{k}(i,:)-Position{k}(1,:))/(Gp-Position{k}(1,:)+1));
            XB(k).temp(i,:)=Position{k}(i,:)+XB(k).dis(i,:).*XB(k).cx(i,:);
            Flag4ub=XB(k).temp(i,:)>ub;
            Flag4lb=XB(k).temp(i,:)<lb;
            XB(k).temp(i,:)=(XB(k).temp(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb; 
            XB(k).tempf(1,i)=func(XB(k).temp(i,:));
            if(XB(k).tempf(1,i)<Fitness{k}(1,i))
                Position{k}(i,:)=XB(k).temp(i,:);
                Fitness{k}(1,i)=XB(k).tempf(1,i);
            end
        end
    end 
    [Position,Fitness]=devide(XB_fitness,XB_position,K,n);% divided into k groups
    [pp,pf,Gp,Gf,P,index,Elite]=compare(Fitness,Position,K,Gp,Gf,pp,pf,n,index,Elite);% update the best individual
    Convergence_curve(t)=Gf;
    t=t+1;
    for k=1:K
        XB_position((k-1)*n+1:n*k,:)=Position{k};
        XB_fitness(1,(k-1)*n+1:n*k)=Fitness{k};
    end
    [Position,Fitness]=devide(XB_fitness,XB_position,K,n);% re-devide into k groups
    %fprintf('BFGO :%05.0fth Function:%02.0f Call:%05.0fth OptimalValue:%04.15f\n',expnum,func_num,1,gbestval(1));
end
end
   
