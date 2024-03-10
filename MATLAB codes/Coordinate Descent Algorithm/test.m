clear all;
clc;

%for i = 1:5    % test 
  %parpool(8)  
    m = 500; 
    n = 500; 
    r = 25; 
    
     %Wt = randn(m,r); 
     %Ht = randn(r,n); 
     %X = max(0,Wt*Ht); 
    %X = full(sprand(m,n,0.2));
    %load CBCL.mat
    load data5001.mat;
    
    [U,S,V] = svds(X,r);
    
    W0 = U;
    H0 = S*V';
    % W0 = randn(m,r); 
    % H0 = randn(r,n); 
    
    
    tic;
    [Wc,Hc,err] = alternatescheme(X,W0,H0,1000); 
    time=toc

%     filename = ['data1000',sprintf('%d',i),'.mat'];
%     save(filename,'X','err','time')
    
    % err : ||X - max(0,WH)||_F / ||X||_F 
    
    semilogy(err); 
    err(end)
%end