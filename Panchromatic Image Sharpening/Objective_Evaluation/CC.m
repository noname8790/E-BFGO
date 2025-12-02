function cc=CC(MS,F) 
d=size(F, 3);
F=double(F);
MS=double(MS);
for i=1:d
    c(i)=corr2(F(:,:,i),MS(:,:,i));
end
cc=(1/3)*(c(1)+c(2)+c(3));
if d>3
cc=(1/4)*(c(1)+c(2)+c(3)+c(4));
end