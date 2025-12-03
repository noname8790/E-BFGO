function f = Costfunc(x,h,u,D)
  alpha=1;lamda=4;gamma=1e4;
  D = D*x';
  s1 = sum((x-h).^2,'all');
  s2 = sum((x-u).^2,'all');
  s3 = sum(D)^2;
  f = alpha*s1 + lamda*s2 + gamma*s3;
end
