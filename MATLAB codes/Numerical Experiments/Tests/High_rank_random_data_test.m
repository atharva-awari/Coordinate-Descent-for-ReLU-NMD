%% Experiment on high-rank synthetic data

clear all;
clc;

   
     m = 200; 
     n = 200; 
     r = 20; 
    
    
    %Generate high rank random matrix
    X = randn(200,200);
    X = max(0,X);
      
    [U,S,V] = svds(X,r);
   
    %Setting parameters for A-NMD, 3B-NMD and CD
    W0 = U;
    H0 = S*V';
    param.Theta0=W0*H0;
    param.W0 = W0;
    param.H0 = H0;
    param.time=1000;
    param.maxit=1000;
    param.stop=500;
    param.tol = 1e-16;
    
%% CD for 1000 iterations
    tic;
    [Wc,Hc,Wc_stop,Hc_stop,err_cd,iter] = alternatescheme(X,W0,H0,param); 
    time_cd=toc;
%% A-NMD for 1000 iterations

    [Theta_anmd,Theta_anmd_500,err_anmd,i_anmd,time_anmd]=A_NMD(X,r,param);
%%  3B-NMD for 1000 iterations
    
     [Theta_3B,err_3B,i_3B,time_3B]=NMD_3B(X,r,param);

%% A-NMD initialized by CD after 500 iterations
    param.Theta0=Wc_stop*Hc_stop;
    param.maxit=500;
    param.time=50;

    [Theta_cd_anmd,Theta_cd_anmd_500,err_cd_anmd,i_cd_anmd,time_cd_anmd]=A_NMD(X,r,param);

 %% CD initialized by A-NMD after 500 iterations
    
    [U_cd,S_cd,V_cd] = svds(Theta_anmd_500,r);
    W0 = U_cd;
    H0 = S_cd*V_cd';

    [Wc_anmd_cd,Hc_anmd_cd,Wc_anmd_cd_stop,Hc_anmd_cd_stop,err_anmd_cd,iter_anmd_cd] = alternatescheme(X,W0,H0,param);




    fprintf('\ntime_CD = %f, iterations_CD = %d, time_ANMD = %f, iterations_ANMD = %d, time_3B = %f, iterations_3B = %d \n', time_cd,iter,time_anmd(end),i_anmd(end),time_3B(end), i_3B(end));

    %% Plotting Graph
    
    figure; 
    set(0, 'DefaultAxesFontSize', 25);
    set(0, 'DefaultLineLineWidth', 2);
    
    semilogy([err_cd(1:param.stop), err_cd_anmd], '--');hold on;
    semilogy([err_anmd(1:param.stop), err_anmd_cd],  ':');
    semilogy(err_cd);
    semilogy(err_anmd);
    semilogy(err_3B);
    
    
    legend('CD(500)+A-NMD(500)', 'A-NMD(500)+CD(500)','CD', 'A-NMD', '3B');
    xlabel('Iterations');
    ylabel('Relative Error');

