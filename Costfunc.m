function cost = Costfunc(x,W_R,W_G,W_B,W_NIR,F_P,F_R,F_G,F_B,F_NIR,I,P,MS_ORG,MS_OBJ,MSWV_US,bandCoeffs)

RESULT(:,:,1) = MSWV_US(:,:,1) + W_R.*(x(1)*F_P+(1-x(1))*F_R).*(P-I);
RESULT(:,:,2) = MSWV_US(:,:,2) + W_G.*(x(2)*F_P+(1-x(2))*F_G).*(P-I);
RESULT(:,:,3) = MSWV_US(:,:,3) + W_B.*(x(3)*F_P+(1-x(3))*F_B).*(P-I);
RESULT(:,:,4) = MSWV_US(:,:,4) + W_NIR.*(x(4)*F_P+(1-x(4))*F_NIR).*(P-I);


    for i=1:size(MS_ORG, 3)
        
    RESULT(:,:,i)=RESULT(:,:,i)*bandCoeffs(i);
    
    end

    cost = 1-QNR(MS_OBJ,P,RESULT,1,1);
    
end


