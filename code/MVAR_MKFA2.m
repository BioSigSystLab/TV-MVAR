%---------------------------------------------------------------------------------------------------------------------
% Copyright (c) 2019 Kyriaki Kostoglou
% PhD Supervisor: Georgios Mitsis, Associate Professor, Bioengineering Department, McGill University, Montreal, Canada
% Biosignals and Systems Analysis Lab
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the % % "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, % distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to % the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.
%---------------------------------------------------------------------------------------------------------------------

%Proposed Kalman filter (improved version) for estimating MVAR models (to be used by the GA). The difference between MVAR_MKFA and MVAR_MKFA2 is that now each model coefficient is assigned a unique forgetting factor for tracking the process noise variance (R1) but also a unique R2 value (noise variance). In MVAR_MKFA, R2 was a diagonal MxM covariance matrix (i.e. each time-series was driven by noise of different variance) but now it is a M*M*pxM*M*p covariance matrix.  

function [J]=MVAR_MKFA2(X,Y,metric,ignore,pmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X:       row vector that contains the hyperparameters optimized by the GA
% Y:       MxN matrix (each time series is placed in a row)
% metric:  fitness function (1 for multivariate AIC, 2 for multivariate BIC)
% ignore:  number of samples to ignore due to initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p=X(1);                             %Model order
[M,N]=size(Y);    
Mp=M*p;    
totPar=M*Mp;                        %Total number of model parameters M*M*p
pin=X(end);                         %P0: Initial value for the covariance matrix P
P = pin*eye(totPar);                %P:  Parameter estimation error covariance matrix
R2=X(2:1+totPar)';                  %R2: Noise variance (however here instead of each input, each parameter is assigned to one value)
R1=zeros(totPar);                   %R1: Initial diagonal process noise covariance matrix (zero)
lam=X(2+M*M*pmax:1+M*M*pmax+totPar);%Forgetting factors assigned to each parameter
LAM=diag(lam);                      %Diagonal matrix containing the forgetting factors assigned to each parameter
yhat=zeros(M,N);                    %Predicted output
th=zeros(totPar,1);                 %M*M*px1 Initial parameter column vector (zero)     
e=zeros(M,N);
phit=zeros(M,totPar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Proposed Kalman Filter (improved)
for k=p+1:N
    temp=Y(:,k-1:-1:k-p);
    temp=temp(:);
    for m=1:M
        phit(m,(m-1)*Mp+1:(m-1)*Mp+Mp)=temp;
    end
    phi=phit';
    yh=phit*th;
    e(:,k)=Y(:,k)-yh;
    epsisq=repmat((e(:,k).^2)',Mp,1);
    R1=LAM.*R1+(eye(totPar)-LAM).*diag(epsisq(:));
    phitp=phit*P;
    rt=phitp*phi;
    rr=repmat(diag(rt)',Mp,1);
    K=(P*phi)./repmat(R2+rr(:),1,M);    
    P=P+R1-K*phitp;
    th=th+K*e(:,k);    
    yhat(:,k)=yh;     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fitness Function
S=cov(e(:,ignore:end)');
ds=det(S);
J=0;
if(ds<0)
    J=inf;
else
    NN=size(Y(:,ignore:end),2);
    if(metric==1)
        J=NN*log(ds)+2*totPar;        %Multivariate AIC
    elseif(metric==2)
        J=NN*log(ds)+log(NN)*totPar;  %Multivariate BIC
    end
end
end



%In order to make a mex file using the matlab coder the inputs type should be as follows:
%X:      double 1xInf
%Y:      double 3x:1000 (note that you can change this based on the size of your data)
%metric: double 1x1
%ignore: double 1x1
%pmax:   double 1x1