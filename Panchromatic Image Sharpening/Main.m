clc;
clear;
close all;
% Loading the images 
addpath WorldView4  %Dataset path
load PAN;
load  MS;
MSWV_db  = double(imgMS);
PANWV_db = double(imgPAN);
MS_ORG   = double(imgMS);
% resizing and upsampling
MSWV_US  = imresize(MSWV_db,  1/4, 'bicubic');
MSWV_US  = imresize(MSWV_US,  4,   'bicubic');
MS_OBJ = MSWV_US;
MSWV_DG  = uint16process(MSWV_US(:,:,1:3));
PANWV_DS = imresize(PANWV_db, 1/4, 'bicubic');
PANWV_US = imresize(PANWV_DS, 4,   'bicubic');
% spectral bands
R   = MSWV_US(:,:,1); 
G   = MSWV_US(:,:,2); 
B   = MSWV_US(:,:,3); 
NIR = MSWV_US(:,:,4); 
% normialization
for i=1:size(MSWV_US,3)
    bandCoeffs(i)      = max(max(MSWV_US(:,:,i)));
    MSWV_US(:,:,i)     = MSWV_US(:,:,i)/bandCoeffs(i);
end
P = PANWV_DS;
panCoeff = max(max(P));
P = P/panCoeff;
% original MS and PAN
figure, imshow(MSWV_DG(:,:,1:3),'border','tight');
figure, imshow(uint16process(PANWV_DS),'border','tight');

% edge detecotr

lamda = 10^-9;
eps   = 10^-10;
F_P   = expEdge(P, lamda,eps);

Red_W = max(max(R));     
R     = R/Red_W;
F_R   = expEdge(R, lamda, eps);

Green_W = max(max(G));    
G       = G/Green_W;
F_G     = expEdge(G, lamda,eps);

Blue_W  = max(max(B));   
B       = B/Blue_W;
F_B     = expEdge(B, lamda,eps);

NIR_W=max(max(NIR));   
NIR=NIR/NIR_W;
F_NIR = expEdge(NIR, lamda,eps);

W = impGradDes(MSWV_US,P);   % Optimal weights for spectral bands

I = W(1)*MSWV_US(:,:,1) + W(2)*MSWV_US(:,:,2) + W(3)*MSWV_US(:,:,3) + W(4)*MSWV_US(:,:,4); 
P = (P-mean(P(:)))*std(I(:))/std(P(:)) + mean(I(:));  % Histogram matching

W_R = 4*R./(R + B + G + NIR);
W_B = 4*B./(R + B + G + NIR);
W_G = 4*G./(R + B + G + NIR);
W_NIR = 4*NIR./(R + B + G + NIR);

addpath Objective_Evaluation
CostFunction=@(x) Costfunc(x,W_R,W_G,W_B,W_NIR,F_P,F_R,F_G,F_B,F_NIR,I,P,MS_ORG,MS_OBJ,MSWV_US,bandCoeffs); % cost function

pop_size = 30; % population size
dim = 4; % problem dimensions
% search boundary
x_min = 0;
x_max = 1;
% max number of iterations
max_it = 50; 

ExpNum = 20;
    
%% optimization algorithm

for expnum=1:ExpNum
[best_pos,bestcost]=PSO_func(CostFunction,dim,pop_size,max_it,x_min,x_max);
best_pos_PSO=best_pos;
% bestcost_PSO=bestcost_PSO+bestcost;bestcost_end_PSO(expnum)=bestcost(end);

[~,best_pos,bestcost]=VPPSO(CostFunction,dim,pop_size,max_it,x_min,x_max);
best_pos_VPPSO=best_pos;
% % bestcost_VPPSO=bestcost_VPPSO+bestcost;bestcost_end_VPPSO(expnum)=bestcost(end);

[~,best_pos,bestcost]=SCA_func(CostFunction,dim,pop_size,max_it,x_min,x_max);
best_pos_SCA=best_pos;
% % bestcost_SCA=bestcost_SCA+bestcost;bestcost_end_SCA(expnum)=bestcost(end);

[best_pos,bestcost]=EDOLSCA(CostFunction,dim,pop_size,max_it,x_min,x_max);
best_pos_EDOLSCA=best_pos;
% % bestcost_EDOLSCA=bestcost_EDOLSCA+bestcost;bestcost_end_EDOLSCA(expnum)=bestcost(end);

[~,best_pos,bestcost]=BFGO(CostFunction,dim,pop_size,max_it,x_min,x_max);
best_pos_BFGO=best_pos;
% % bestcost_BFGO=bestcost_BFGO+bestcost;bestcost_end_BFGO(expnum)=bestcost(end);

[~,best_pos,bestcost]=EBFGO(CostFunction,dim,pop_size,max_it,x_min,x_max);
best_pos_EBFGO=best_pos;
% bestcost_EBFGO=bestcost_EBFGO+bestcost;bestcost_end_EBFGO(expnum)=bestcost(end);
end

% PCA method
RESULT_PCA = PCA(MS_ORG, PANWV_DS);

% plot(1:MaxIt,bestcost_PSO,'b','LineWidth',1.2);hold on;
% plot(1:MaxIt,bestcost_VPPSO,'color',[0.5,0.5,0.5],'LineWidth',1.2);hold on;
% plot(1:MaxIt,bestcost_SCA,'g','LineWidth',1.2);hold on;
% plot(1:MaxIt,bestcost_EDOLSCA,'y','LineWidth',1.2);hold on;
% plot(1:MaxIt,bestcost_BFGO,'color',[1,0.75,0.8],'LineWidth',1.2);hold on;
% plot(1:MaxIt,bestcost_EBFGO,'r','LineWidth',1.2);hold off;
% set(gca,'Fontname','Times New Roman','FontSize',10);
% legend('PSO','VPPSO','SCA','EDOLSCA','BFGO','EBFGO');
% title('CostFunction Value For Each Function');
% xlabel('Iteration');
% ylabel('CostFunction Value');

% disp(['Best Cost of PSO = ' num2str(min(bestcost_end_PSO))]);
% disp(['Best Cost of VPPSO = ' num2str(min(bestcost_end_VPPSO))]);
% disp(['Best Cost of SCA = ' num2str(min(bestcost_end_SCA))]);
% disp(['Best Cost of EDOLSCA = ' num2str(min(bestcost_end_EDOLSCA))]);
% disp(['Best Cost of BFGO = ' num2str(min(bestcost_end_BFGO))]);
% disp(['Best Cost of EBFGO = ' num2str(min(bestcost_end_EBFGO))]);


%% calculate fusion result

RESULT_PSO(:,:,1) = MSWV_US(:,:,1) + W_R.*(best_pos_PSO(1)*F_P + (1-best_pos_PSO(1))*F_R).*(P-I);
RESULT_PSO(:,:,2) = MSWV_US(:,:,2) + W_G.*(best_pos_PSO(2)*F_P + (1-best_pos_PSO(2))*F_G).*(P-I);
RESULT_PSO(:,:,3) = MSWV_US(:,:,3) + W_B.*(best_pos_PSO(3)*F_P + (1-best_pos_PSO(3))*F_B).*(P-I);
RESULT_PSO(:,:,4) = MSWV_US(:,:,4) + W_NIR.*(best_pos_PSO(4)*F_P + (1-best_pos_PSO(4))*F_NIR).*(P-I);

RESULT_VPPSO(:,:,1) = MSWV_US(:,:,1) + W_R.*(best_pos_VPPSO(1)*F_P + (1-best_pos_VPPSO(1))*F_R).*(P-I);
RESULT_VPPSO(:,:,2) = MSWV_US(:,:,2) + W_G.*(best_pos_VPPSO(2)*F_P + (1-best_pos_VPPSO(2))*F_G).*(P-I);
RESULT_VPPSO(:,:,3) = MSWV_US(:,:,3) + W_B.*(best_pos_VPPSO(3)*F_P + (1-best_pos_VPPSO(3))*F_B).*(P-I);
RESULT_VPPSO(:,:,4) = MSWV_US(:,:,4) + W_NIR.*(best_pos_VPPSO(4)*F_P + (1-best_pos_VPPSO(4))*F_NIR).*(P-I);

RESULT_SCA(:,:,1) = MSWV_US(:,:,1) + W_R.*(best_pos_SCA(1)*F_P + (1-best_pos_SCA(1))*F_R).*(P-I);
RESULT_SCA(:,:,2) = MSWV_US(:,:,2) + W_G.*(best_pos_SCA(2)*F_P + (1-best_pos_SCA(2))*F_G).*(P-I);
RESULT_SCA(:,:,3) = MSWV_US(:,:,3) + W_B.*(best_pos_SCA(3)*F_P + (1-best_pos_SCA(3))*F_B).*(P-I);
RESULT_SCA(:,:,4) = MSWV_US(:,:,4) + W_NIR.*(best_pos_SCA(4)*F_P + (1-best_pos_SCA(4))*F_NIR).*(P-I);

RESULT_EDOLSCA(:,:,1) = MSWV_US(:,:,1) + W_R.*(best_pos_EDOLSCA(1)*F_P + (1-best_pos_EDOLSCA(1))*F_R).*(P-I);
RESULT_EDOLSCA(:,:,2) = MSWV_US(:,:,2) + W_G.*(best_pos_EDOLSCA(2)*F_P + (1-best_pos_EDOLSCA(2))*F_G).*(P-I);
RESULT_EDOLSCA(:,:,3) = MSWV_US(:,:,3) + W_B.*(best_pos_EDOLSCA(3)*F_P + (1-best_pos_EDOLSCA(3))*F_B).*(P-I);
RESULT_EDOLSCA(:,:,4) = MSWV_US(:,:,4) + W_NIR.*(best_pos_EDOLSCA(4)*F_P + (1-best_pos_EDOLSCA(4))*F_NIR).*(P-I);

RESULT_BFGO(:,:,1) = MSWV_US(:,:,1) + W_R.*(best_pos_BFGO(1)*F_P + (1-best_pos_BFGO(1))*F_R).*(P-I);
RESULT_BFGO(:,:,2) = MSWV_US(:,:,2) + W_G.*(best_pos_BFGO(2)*F_P + (1-best_pos_BFGO(2))*F_G).*(P-I);
RESULT_BFGO(:,:,3) = MSWV_US(:,:,3) + W_B.*(best_pos_BFGO(3)*F_P + (1-best_pos_BFGO(3))*F_B).*(P-I);
RESULT_BFGO(:,:,4) = MSWV_US(:,:,4) + W_NIR.*(best_pos_BFGO(4)*F_P + (1-best_pos_BFGO(4))*F_NIR).*(P-I);

RESULT_EBFGO(:,:,1) = MSWV_US(:,:,1) + W_R.*(best_pos_EBFGO(1)*F_P + (1-best_pos_EBFGO(1))*F_R).*(P-I);
RESULT_EBFGO(:,:,2) = MSWV_US(:,:,2) + W_G.*(best_pos_EBFGO(2)*F_P + (1-best_pos_EBFGO(2))*F_G).*(P-I);
RESULT_EBFGO(:,:,3) = MSWV_US(:,:,3) + W_B.*(best_pos_EBFGO(3)*F_P + (1-best_pos_EBFGO(3))*F_B).*(P-I);
RESULT_EBFGO(:,:,4) = MSWV_US(:,:,4) + W_NIR.*(best_pos_EBFGO(4)*F_P + (1-best_pos_EBFGO(4))*F_NIR).*(P-I);

for i = 1:size(MS_ORG, 3)
    RESULT_PSO(:,:,i)= RESULT_PSO(:,:,i)*bandCoeffs(i);
    RESULT_VPPSO(:,:,i)= RESULT_VPPSO(:,:,i)*bandCoeffs(i);
    RESULT_SCA(:,:,i)= RESULT_SCA(:,:,i)*bandCoeffs(i);
    RESULT_EDOLSCA(:,:,i)= RESULT_EDOLSCA(:,:,i)*bandCoeffs(i);
    RESULT_BFGO(:,:,i)= RESULT_BFGO(:,:,i)*bandCoeffs(i);
    RESULT_EBFGO(:,:,i)= RESULT_EBFGO(:,:,i)*bandCoeffs(i);
end

% display the fusion results
figure, imshow(uint16process(RESULT_PCA(:,:,1:3)),'border','tight');
figure, imshow(uint16process(RESULT_PSO(:,:,1:3)),'border','tight');
figure, imshow(uint16process(RESULT_VPPSO(:,:,1:3)),'border','tight');
figure, imshow(uint16process(RESULT_SCA(:,:,1:3)),'border','tight');
figure, imshow(uint16process(RESULT_EDOLSCA(:,:,1:3)),'border','tight');
figure, imshow(uint16process(RESULT_BFGO(:,:,1:3)),'border','tight');
figure, imshow(uint16process(RESULT_EBFGO(:,:,1:3)),'border','tight');
% evaluating results

DS_PCA      = D_s(MS_OBJ,P,RESULT_PCA);
Dlambda_PCA = D_lambda(MS_OBJ,RESULT_PCA);
QNR_PCA     = QNR(MS_OBJ,P,RESULT_PCA,1,1);
CC_PCA      = CC(MS_OBJ,RESULT_PCA);
T_PCA       = table(DS_PCA, Dlambda_PCA, QNR_PCA, CC_PCA ...
    ,'VariableNames', {'D_s', 'D_lambda', 'QNR', 'CC'});


DS_PSO      = D_s(MS_OBJ,P,RESULT_PSO);
Dlambda_PSO = D_lambda(MS_OBJ,RESULT_PSO);
QNR_PSO     = QNR(MS_OBJ,P,RESULT_PSO,1,1);
CC_PSO      = CC(MS_OBJ,RESULT_PSO);
T_PSO       = table(DS_PSO, Dlambda_PSO, QNR_PSO, CC_PSO ...
    ,'VariableNames', {'D_s', 'D_lambda', 'QNR', 'CC'});


DS_VPPSO      = D_s(MS_OBJ,P,RESULT_VPPSO);
Dlambda_VPPSO = D_lambda(MS_OBJ,RESULT_VPPSO);
QNR_VPPSO     = QNR(MS_OBJ,P,RESULT_VPPSO,1,1);
CC_VPPSO      = CC(MS_OBJ,RESULT_VPPSO);
T_VPPSO       = table(DS_VPPSO, Dlambda_VPPSO, QNR_VPPSO, CC_VPPSO ...
    ,'VariableNames', {'D_s', 'D_lambda', 'QNR', 'CC'});


DS_SCA      = D_s(MS_OBJ,P,RESULT_SCA);
Dlambda_SCA = D_lambda(MS_OBJ,RESULT_SCA);
QNR_SCA     = QNR(MS_OBJ,P,RESULT_SCA,1,1);
CC_SCA      = CC(MS_OBJ,RESULT_SCA);
T_SCA       = table(DS_SCA, Dlambda_SCA, QNR_SCA, CC_SCA ...
    ,'VariableNames', {'D_s', 'D_lambda', 'QNR', 'CC'});


DS_EDOLSCA      = D_s(MS_OBJ,P,RESULT_EDOLSCA);
Dlambda_EDOLSCA = D_lambda(MS_OBJ,RESULT_EDOLSCA);
QNR_EDOLSCA     = QNR(MS_OBJ,P,RESULT_EDOLSCA,1,1);
CC_EDOLSCA      = CC(MS_OBJ,RESULT_EDOLSCA);
T_EDOLSCA       = table(DS_EDOLSCA, Dlambda_EDOLSCA, QNR_EDOLSCA, CC_EDOLSCA ...
    ,'VariableNames', {'D_s', 'D_lambda', 'QNR', 'CC'});


DS_BFGO      = D_s(MS_OBJ,P,RESULT_BFGO);
Dlambda_BFGO = D_lambda(MS_OBJ,RESULT_BFGO);
QNR_BFGO     = QNR(MS_OBJ,P,RESULT_BFGO,1,1);
CC_BFGO      = CC(MS_OBJ,RESULT_BFGO);
T_BFGO       = table(DS_BFGO, Dlambda_BFGO, QNR_BFGO, CC_BFGO ...
    ,'VariableNames', {'D_s', 'D_lambda', 'QNR', 'CC'});


DS_EBFGO      = D_s(MS_OBJ,P,RESULT_EBFGO);
Dlambda_EBFGO = D_lambda(MS_OBJ,RESULT_EBFGO);
QNR_EBFGO     = QNR(MS_OBJ,P,RESULT_EBFGO,1,1);
CC_EBFGO      = CC(MS_OBJ,RESULT_EBFGO);
T_EBFGO       = table(DS_EBFGO, Dlambda_EBFGO, QNR_EBFGO, CC_EBFGO ...
    ,'VariableNames', {'D_s', 'D_lambda', 'QNR', 'CC'});

% RESULT=[T_EBFGO;T_BFGO;T_VPPSO;T_PSO;T_EDOLSCA;T_SCA];
% algorithm_names = {'EBFGO', 'BFGO','VPPSO', 'PSO', 'EDOLSCA','SCA'  };
% RESULT.Properties.RowNames = algorithm_names;
% writetable(RESULT,'result.xlsx');

