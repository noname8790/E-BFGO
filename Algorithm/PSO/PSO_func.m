function [gbest,gbestval]= PSO_func(func,dim,popsize,maxgen,xmin,xmax)
xmin=xmin.*ones(1,dim); 
xmax=xmax.*ones(1,dim);
c1=2.05;
c2=2.05;
w=0.9;
alpha=1;
gbestval=zeros(1,maxgen);
vmax=0.5.*(xmax-xmin).*ones(popsize,dim);
vmin=-vmax;
v=vmin+2.*vmax.*rand(popsize,dim);%initialize the velocity of the particles
x=xmin+(xmax-xmin).*rand(popsize,dim);
pbest=x;
for i=1:popsize
pbestval(i)=func(x(i,:));%initialize the Pbest and the Pbest's fitness value
end
[gbestval(1),gbestid]=min(pbestval);
gbest=pbest(gbestid,:);%initialize the gbest and the gbest's fitness value
for t=2:maxgen
    gbestrep=repmat(gbest,popsize,1);%update the gbest
    v=w*v+c1.*rand(popsize,dim).*(pbest-x)+c2.*rand(popsize,dim).*(gbestrep-x);
    v=(v>vmax).*vmax+(v<=vmax).*v;
    v=(v<vmin).*vmin+(v>=vmin).*v;
    x=x+alpha*v;
    x=((x>=xmin)&(x<=xmax)).*x+(x<xmin).*(xmin+0.25.*(xmax-xmin).*rand(popsize,dim))+(x>xmax).*(xmax-0.25.*(xmax-xmin).*rand(popsize,dim));    
    for i=1:popsize
    e(i)=func(x(i,:));
    end
    tmp=(pbestval<e);
    temp=repmat(tmp',1,dim);
    pbest=temp.*pbest+(1-temp).*x;
    pbestval=tmp.*pbestval+(1-tmp).*e;%update the Pbest
    [gbestval(t),gbestid]=min(pbestval);
    gbest=pbest(gbestid,:);    
end
end


