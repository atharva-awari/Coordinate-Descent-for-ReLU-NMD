%% Experiment on Synthetic data

clear all;
clc;

  m = 50;n=50; r=8;
   
   %Generate a synthetic matrix
   Z = randn(m,n);
   X = max(0,Z);
      
    [U,S,V] = svds(X,r);
    
    % Parameter Setting for CD and A-NMD
    W0 = U;
    H0 = S*V';
    param.Theta0=W0*H0;
    param.W0 = W0;
    param.H0 = H0;
    param.time=1000;
    param.maxit=200;
    param.tol = 1e-4;
    param.tolerr=1e-4;
    
%%  CD for 1000 iterations
    tic;
    [Wc,Hc,Wc_stop,Hc_stop,err_cd,iter] = alternatescheme(X,W0,H0,param); 
    time_cd=toc;
%%  A-NMD for 1000 iterations

    [Theta_anmd,Theta_anmd_stop,err_anmd,i_anmd,time_anmd]=A_NMD(X,r,param);
    
    [Theta_3B,err_3B,i_3B,time_3B]=NMD_3B(X,r,param);


    %% Plotting the graph
    
    figure; 
    set(0, 'DefaultAxesFontSize', 25);
    set(0, 'DefaultLineLineWidth', 3);
    

    semilogy(err_cd); hold on;
    semilogy(err_anmd,':');
    semilogy(err_3B,'-.');
    
    
    legend('CD', 'A-NMD', '3B');
    xlabel('Iterations');
    ylabel('Relative Error');