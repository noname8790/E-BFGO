function [Gf,Gp,Convergence_curve] = EBFGO(func,D,N,T,lb,ub)
%% initialize bamboo population
rand('twister',mod(floor(now*8640000),2^31-1));
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
XB_lifetime=ones(1,N);% initialize the age of the individual
alpha=0.2; % elimination rate
base_nutrient = 2;% initial nutrient index
stalk_time=3;% the age when bamboo shoots emerge
dead_time=8;% life of bamboo
for k=1:K
    pf{k}=Inf;
end
index=1;
[Position,Fitness,lifetime]=devide(XB_fitness,XB_position,XB_lifetime,K,n);% devided into k groups
[pp,pf,Gp,Gf,P,index,Elite]=initial_compare(Fitness,Position,K,Gp,Gf,pp,pf,n,index);% update the best individual
Convergence_curve(1)=Gf;
t=2;
BW_Position=Position;% bamboo whip
BW_Fitness=Fitness;
BS_Position=Position;% bamboo shoot
BS_Fitness=Fitness;
%% iteration
while t<T+1
    if(t>10&&Convergence_curve(t-1)==Convergence_curve(t-3))
        [good_position,poor_position]=fliped(BW_Position,BW_Fitness,K,Gp,Gf,pp,pf,n,Elite);
        for k=1:K
            for i=1:n-2*n/5
                good_fitness{k}(1,i)=func(good_position{k}(i,:));
                if(good_fitness{k}(1,i)<BW_Fitness{k}(1,i))
                    BW_Position{k}(i,:)=good_position{k}(i,:);
                    BW_Fitness{k}(1,i)=good_fitness{k}(1,i);
                end
            end
            for i=1:2*n/5
                poor_fitness{k}(1,i)=func(poor_position{k}(i,:));
            end
            BW_Position{k}(n-2*n/5+1:n,:)=poor_position{k};
            BW_Fitness{k}(1,n-2*n/5+1:n)=poor_fitness{k};
        end
        [pp,pf,Gp,Gf,P,index,Elite]=compare(BW_Fitness,BW_Position,K,Gp,Gf,pp,pf,n,index,Elite);% update the best individual
    end   
    if(sum(P)==K)
        for k=1:K
        BW_Position{k}=XB_position((k-1)*n+1:k*n,:); 
        BW_Fitness{k}=XB_fitness(1,(k-1)*n+1:k*n); 
        lifetime{k}=XB_lifetime(1,(k-1)*n+1:k*n);
        end  
    end

%% extended bamboo whip extension stage
    a=1.2-(t/T);
    soil_nutrient = base_nutrient * exp(-0.5*(t/T)); % nutrient index decay
    for k=1:K
        C{k}=mean(BW_Position{k},1);
        for i=1:n
            r1=rand();
            % adjust weights based on soil nutrients and individual age
            extend_factor=(0.7*soil_nutrient+0.3*(stalk_time/lifetime{k}(1,i)))^a;
            if(r1<0.4)
               cossita=dot(Gp,BW_Position{k}(i,:))/norm1(Gp)/norm1(BW_Position{k}(i,:));
               BW_Position{k}(i,:)=Gp+a*extend_factor*(Gp-BW_Position{k}(i,:))*cossita;
            elseif(r1>0.4&&r1<0.7) 
               cossita=dot(pp{k},BW_Position{k}(i,:))/norm1(pp{k})/norm1(BW_Position{k}(i,:));
               BW_Position{k}(i,:)=pp{k}+a*extend_factor*(pp{k}-BW_Position{k}(i,:))*cossita;
            else
               cossita=dot(C{k},BW_Position{k}(i,:))/norm1(C{k})/norm1(BW_Position{k}(i,:));
               BW_Position{k}(i,:)=C{k}+a*extend_factor*(C{k}-BW_Position{k}(i,:))*cossita;
            end
            BW_Position{k}(i,:)=checkbound(BW_Position{k}(i,:),lb,ub);
        end  
        % one local search
        for i=1:n 
            perturbation=unifrnd(lb,ub,1,D);
            new_position=BW_Position{k}(i,:) + perturbation/((ub-lb)*5);
            new_position=checkbound(new_position,lb,ub);
            new_fitness=func(new_position);
            if (new_fitness<BW_Fitness{k}(1,i))
                BW_Position{k}(i,:)=new_position;
                BW_Fitness{k}(1,i)=new_fitness;
            end
        end
    end
    [pp,pf,Gp,Gf,P,index,Elite]=compare(BW_Fitness,BW_Position,K,Gp,Gf,pp,pf,n,index,Elite);% update the best individual

%% extended bamboo shoot growth stage
    % growth level ranking
    [~,rank_c] = sort(cellfun(@mean, BS_Fitness), 'descend');
    sorted_c(rank_c) = 1:length(rank_c);
    % relate growth model parameters
    H_max = ub; k_growth = 1.5; 
    branch_angle = 34.14*pi/180; % the angle between bamboo branches and stalks (golden ratio point)
    rot_matrix = [cos(branch_angle) -sin(branch_angle); % branches rotation matrix
                 sin(branch_angle)  cos(branch_angle)]; 
    global_light=randn()*2;
    for k = 1:K
        light_factor = global_light *sorted_c(k)/K; % lighting available to each cluster
        for i=1:n
            if lifetime{k}(1,i)>=stalk_time
                % bamboo stalk growth length[logistic]
                growth_main = k_growth*(1-norm(BS_Position{k}(i,:)/H_max,2)).* ((Gp-BS_Position{k}(i,:))/H_max).*( light_factor+soil_nutrient);
                grown_pos = BS_Position{k}(i,:) + growth_main;
                grown_pos = checkbound(grown_pos,lb,ub);
                random_dim = randperm(D); 
                % bamboo branch norm length
                branch_norm = norm(growth_main/max(growth_main),2) * (1-(lifetime{k}(1,i)+1-stalk_time)/(dead_time+1-stalk_time));% growth decay
                grown_pos_with_branch = grown_pos;
                j=2;
                % grow branches in random directions
                while j<=D
                r2=rand();
                if r2<=0.3
                grown_pos_with_branch(random_dim(j-1:j)) = grown_pos_with_branch(random_dim(j-1:j)) + (rot_matrix * [abs(randn()*branch_norm); 0])';
                elseif r2>0.3&&r2<=0.6
                grown_pos_with_branch(random_dim(j-1:j)) = grown_pos_with_branch(random_dim(j-1:j)) + (rot_matrix * [0; abs(randn()*branch_norm)])';
                end
                j = j+2;
                end
                grown_pos_with_branch = checkbound(grown_pos_with_branch,lb,ub);
                r3=rand();
                % bamboo stalk without branches
                if r3<=0.8
                fit = func(grown_pos_with_branch);
                if fit<BS_Fitness{k}(1,i)
                BS_Position{k}(i,:) = grown_pos_with_branch;
                BS_Fitness{k}(1,i) = fit;
                end
                else
                fit = func(grown_pos);
                if fit<BS_Fitness{k}(1,i)
                BS_Position{k}(i,:) = grown_pos;
                BS_Fitness{k}(1,i) = fit;
                end
                end
            end
        end
    end
%% eliminate and renew the population
    [BW_Position,BW_Fitness,BS_Position,BS_Fitness,lifetime] = elimitate(BW_Position,BW_Fitness,BS_Position,BS_Fitness,lifetime,stalk_time,dead_time, ...
        func,n,D,K,lb,ub,alpha);
    [pp,pf,Gp,Gf,P,index,Elite] = compare(BS_Fitness,BS_Position,K,Gp,Gf,pp,pf,n,index,Elite);% update the best individual
    Convergence_curve(t)=Gf;
    t=t+1;
    for k=1:K
        lifetime{k}=lifetime{k}+1;
        XB_lifetime(1,(k-1)*n+1:n*k)=lifetime{k};
        XB_position((k-1)*n+1:n*k,:)=BW_Position{k};
        XB_fitness(1,(k-1)*n+1:n*k)=BW_Fitness{k};
    end
    [BW_Position,BW_Fitness,lifetime]=devide(XB_fitness,XB_position,XB_lifetime,K,n);% re-devide into k groups
end
end

function [Position] = checkbound(Position,lb,ub)
    % transboundary reflect
    over_ub = Position > ub;
    over_lb = Position < lb;
    Position(over_ub) = 2*ub - Position(over_ub);
    Position(over_lb) = 2*lb - Position(over_lb);
end

function [BW_Position,BW_Fitness,BS_Position,BS_Fitness,lifetime] = elimitate(BW_Position,BW_Fitness,BS_Position,BS_Fitness, ...
    lifetime,stalk_time,dead_time,func,n,D,K,lb,ub,alpha)
    % determine the number of eliminations by eliminate rate
    el_num = round(alpha * n); 
    % logistic mapping-initialization
    init_sum =K*n;
    Z = zeros(init_sum, D);
    Z(1,:) = rand(1, D);
    for i=2:init_sum
    Z(i,:) = 4.0*Z(i-1,:).*(1-Z(i-1,:));
    end
    i=1;
    for k=1:K
        % eliminate worse individuals
        [~, rank] = sort(BS_Fitness{k}, 'descend'); 
        el_idxs = rank(1:el_num);
        for j = 1:length(el_idxs)
            el_idx = el_idxs(j);
            BW_Position{k}(el_idx, :) = Z(i,:) * (ub - lb) + lb; 
            BW_Position{k}(el_idx, :) = checkbound(BW_Position{k}(el_idx, :),lb,ub);
            BW_Fitness{k}(1, el_idx) = func(BW_Position{k}(el_idx, :));
            if lifetime{k}(1,j)>=stalk_time
            BS_Position=BW_Position;BS_Fitness=BW_Fitness;
            end
            lifetime{k}(1,el_idx) = 1;
            i = i+1;
        end
        % eliminate older individuals
        for j=1:n
            if lifetime{k}(1,j)>=dead_time
            BW_Position{k}(j,:) = Z(i,:) * (ub - lb) + lb; 
            BW_Position{k}(j,:) = checkbound(BW_Position{k}(j,:),lb,ub);
            BW_Fitness{k}(1, j) = func(BW_Position{k}(j, :));
            BS_Position=BW_Position;BS_Fitness=BW_Fitness;
            lifetime{k}(1,j) = 1;
            i=i+1;
            end    
        end
    end
end
