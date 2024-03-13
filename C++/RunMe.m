clear all; close all; clc;

%To run this code, we have to compile ACD_NMD.cpp
% In Octave :
% - write this command :
%    mkoctfile --mex ACD_NMD.cpp -lgsl -lgslcblas
%
%In Matlab :
% - download gsl library  : https://www.gnu.org/software/gsl/
%    (only compatible with Linux)
% - write this command:
%    mex -I/path/to/gsl/include -L/path/to/gsl/lib -lgsl -lgslcblas Alt_NMD.cpp


%Fix the size of the matrix X
m=100;n=100;

%Choose the approximation rank
r=32;

%Generate the synthetic matrix
Wt = randn(m,r);
Ht = randn(r,n);
X = max(0,Wt*Ht);


%for i=1:5
%  filename = ['data1000',sprintf('%d',i),'.mat'];
%  load(filename);

%[U,S,V] = svds(X,r);

%Run the codes with the default choice for the parameters
%Or Change the parameters below
%param.W0=  U;
%param.H0=  S*V';
%param.maxiter = ;
#param.tol =  ;
%param.beta = ;
%param.eta = ;
%param.gamma = ;
%param.gamma0 = ;

%ACD_NMD
tic;
[Theta,err_cpp,t] = alternatescheme_acc(X,r);
time = toc

%Resize the vector err and time
pos  = abs(err_cpp) > 1e-16;
err_cpp = err_cpp(pos);
t = t(pos);
 iterations = size(err_cpp)(1)


%%Plot the error per iteration
figure
title('Alternate scheme for NMD accelerated');
subplot(1,2,1);
semilogy(err_cpp, 'LineWidth',1.5);
grid on;
xlabel('Iterations')


%%Plot the error as function of time
subplot(1,2,2);
semilogy(t, err_cpp, 'LineWidth',1.5);
grid on;
xlabel('Time (s.)')

%end

%%Compare errors from MATLAB and Octave
%load err1_sparse.mat;
%error_matlab = err;

% Specify the tolerance (e.g., 1e-5 for 5 digits of accuracy)
%tolerance = 1e-5;

% Compare arrays with specified tolerance
%result = all(abs(error_matlab - err_cpp) < tolerance);

% Display the result
%if all(result)==1
%    disp('Arrays are approximately equal within 5 digits of accuracy.');
%else
%    disp('Arrays are not approximately equal within 5 digits of accuracy.');
%end

