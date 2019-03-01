%---------------------------------------------------------------------------------------------------------------------
% Copyright (c) 2019 Kyriaki Kostoglou
% PhD Supervisor: Georgios Mitsis, Associate Professor, Bioengineering Department, McGill University, Montreal, Canada
% Biosignals and Systems Analysis Lab
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the % % "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, % distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to % the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.
%---------------------------------------------------------------------------------------------------------------------

%Conventional Kalman Filter for estimating MVAR models (to be used by the GA)

function [J]=MVAR_KF(X,Y,metric,ignore)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X:       row vector that contains the hyperparameters optimized by the GA
% Y:       MxN matrix (each time series is placed in a row)
% metric:  fitness function (1 for multivariate AIC, 2 for multivariate BIC)
% ignore:  number of samples to ignore due to initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p=X(1);             %Model order
[M,N]=size(Y);    
Mp=M*p;    
totPar=M*Mp;        %Total number of model parameters M*M*p
pin=X(end);         %P0: Initial value for the covariance matrix P
P = pin*eye(Mp);    %P:  Parameter estimation error covariance matrix
R2=X(2);            %R2: Measurement noise variance
R1=X(3)*eye(Mp);    %R1: Diagonal process noise covariance matrix
yhat=zeros(M,N);    %Predicted output
th=zeros(M,Mp);     %MxMp Initial parameter vector (zero)     
e=zeros(M,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kalman Filter (see Eqs.1-7)
for k=p+1:N
    temp= Y(:,k-1:-1:k-p);
    phi=temp(:);
    phit=phi';
    yh=th*phi;
    e(:,k)=Y(:,k)-yh;
    phitP=phit*P;
    rt=phitP*phi;
    K=(P*phi)/(R2+rt); 
    P=P+R1-K*phitP;
    th=th+(K*e(:,k)')';  
    yhat(:,k)=yh;     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fitness Function
J=0;
S=cov(e(:,ignore:end)');
ds=det(S);
if(ds<0)
    J=inf;
else
    NN=size(Y(:,ignore:end),2);
    if(metric==1)
        J=NN*log(ds)+2*totPar;       %Multivariate AIC
    elseif(metric==2)
        J=NN*log(ds)+log(NN)*totPar; %Multivariate BIC
    end
end





%In order to make a mex file using the matlab coder the inputs type should be as follows:
%X:      double 1xInf
%Y:      double 3x:1000 (note that you can change this based on the size of your data)
%metric: double 1x1
%ignore: double 1x1
