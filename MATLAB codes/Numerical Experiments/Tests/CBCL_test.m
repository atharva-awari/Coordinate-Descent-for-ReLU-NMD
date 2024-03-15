%% Experiment on sparse images obtained through NMF

clear all
close all
clc;

% Add paths
cd('../'); 
Install; 

%Parameters setting
param.maxit=10000; param.tol=1.e-16; param.tolerr = 0; param.time=70; param.stop=5000;

%A-NMD acceleration parameters
param.beta=0.7; param.eta=0.4; param.gamma=1.1; param.gamma_bar=1.05;

%Load the CBCL 19x19 faces data set
load CBCL 
s = 100;      %approximation rank of the NMF factors
%Compute NMF using minvol-NMF
[V,U] = minvolNMF(X',s); 
U = U'; [m,n] = size(U);

%Compute the best V
V = NNLS(U,X); 

%NMF error evaluation 
err_NMF_real=norm(X-U*V,'fro')/norm(X,'fro');

r=20;
[U_svd,S_svd,V_svd] = svds(U,r);
    
    W0 = U_svd;
    H0 = S_svd*V_svd';
    param.Theta0=W0*H0;
      
    param.time=15000;
    param.maxit=10000;
    param.stop=5000;
    param.tol = 1e-16;
%% CD for 10000 iterations
  tic;
    [Wc,Hc,Wc_stop,Hc_stop,err_cd,iter] = alternatescheme(U,W0,H0,param); 
  time_cd=toc;

%% A-NMD for 10000 iterations
[U_ANMD_neg,U_ANMD_neg_stop, err_acc_and,it_acc_and,t_and_acc]=A_NMD(U,r,param);
err_acc_and(end)
U_ANMD=max(0,U_ANMD_neg);
time_ANMD=t_and_acc(end)/it_acc_and;
err_ANMD=norm(U-U_ANMD,'fro')/norm(U,'fro');            %NMD error
V_ANMD = NNLS(U_ANMD,X);                                   %Best V
err_ANMD_nmf=norm(X-U_ANMD*V_ANMD,'fro')/norm(X,'fro'); %NMF error



%% CD initialized by A-NMD after 5000 iterations

[P,Q,R] = svds(U_ANMD_neg_stop, r);
W0 = P;
H0 = Q*R';
param.maxit=50;

[W_anmd_cd,H_anmd_cd,Wc_anmd_cd_stop,Hc_anmd_cd_stop,err_anmd_cd,iter_anmd_cd] = alternatescheme(U,W0,H0,param);

%% A-NMD initialized by CD after 5000 iterations

    param.Theta0=Wc_stop*Hc_stop;
    param.maxit=50;
    param.time=50;

    [U_ANMD_neg_cd_anmd,U_ANMD_neg_stop_cd_anmd, err_acc_and_cd_anmd,it_acc_and_cd_anmd,t_and_acc_cd_anmd]=A_NMD(U,r,param);

%% Figure

figure; 
    set(0, 'DefaultAxesFontSize', 25);
    set(0, 'DefaultLineLineWidth', 2);
    
    semilogy([err_cd(1:param.stop), err_acc_and_cd_anmd], '--');hold on;
    semilogy([err_acc_and(1:param.stop), err_anmd_cd],  ':');
    semilogy(err_cd);
    semilogy(err_acc_and);
    
   
    legend('CD+A-NMD', 'A-NMD+CD','CD', 'A-NMD');
    xlabel('Iterations');
    ylabel('Relative Error');