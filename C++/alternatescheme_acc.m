function [Theta, err, time] = alternatescheme_acc(X,r,param)

%Compute an approximate solution
%of the following non linear matrix decomposition problem (NMD)
%
%    min_{X,W,H} ||X - max(0,WH)||_F  s.t. rank(WH)<=r
%
%This algorithm is an alternate scheme for NMD using Coordinate descent and extrapolation for the acceleration
%and compiling in C++
%
%****** Input ******
%   X       : m-by-n matrix, sparse and non negative
%   param   : structure, containing the parameter of the model
%       .W0 = initialization of the variable W (default: randn)
%       .H0 = initialization of the variable H (default: randn)
%       .maxiter  = maximum number of iterations (default: 100)
%       .tol    = tolerance on the relative error (default: 1.e-4)
%       Coefficients of extrapolation
%       .beta   (default: 0.25)
%       .eta    (default: 1.5)
%       .gamma  (default: 1.05)
%       .gamma0 (default: 1.005)
%
% ****** Output ******
%   Theta   : m-by-n matrix, approximate solution of
%             min_{Theta}||X-max(0,Theta)||_F^2  s.t. rank(Theta)=r
%   err     : vector containing evolution of relative error along
%             iterations ||X-max(0,Theta)||_F / || X ||_F
%   time    : vector containing time counter along iterations

    [m,n]=size(X);
    if nargin < 3
      param = [];
    end

    if ~isfield(param,'W0') || ~isfield(param,'H0')
      param.W0=randn(m,r); param.H0=randn(r,n);
    end

    if ~isfield(param,'maxiter')
      param.maxiter = 500;
    end
    if ~isfield(param,'tol')
      param.tol = 1.e-4;
    end
    if ~isfield(param,'beta')
      param.beta = 0.25;
    end
    if ~isfield(param,'eta')
      param.eta = 1.5;
    end
    if ~isfield(param,'gamma')
      param.gamma= 1.05;
    end
    if ~isfield(param,'gamma0')
      param.gamma0 = 1.005;
    end

    time = zeros(param.maxiter + 1,1);

    % Scale W,H
    WH = param.W0*param.H0;
    Z = max(0,WH);
    optscal = sum(sum(Z.*X))/sum(Z(:).^2);
    param.W0 = optscal*param.W0;
    WH = optscal*WH;

    %Compute the error
    normX = norm(X,'fro');
    err = zeros(param.maxiter + 1,1);
    e = [norm(X-optscal*Z,'fro')];
    e = e/normX;

    %Call the mexFile
%     mex ACD_NMD.cpp;
    [H, W, err, t] = ACD_NMD(X, param.W0, param.H0, param.maxiter, WH, normX, param.beta, param.eta, param.gamma, param.gamma0, e, param.tol);

    time(2:param.maxiter+1) = t;

    Theta = max(0,W*H);
end
