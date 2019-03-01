%---------------------------------------------------------------------------------------------------------------------
% Copyright (c) 2019 Kyriaki Kostoglou
% PhD Supervisor: Georgios Mitsis, Associate Professor, Bioengineering Department, McGill University, Montreal, Canada
% Biosignals and Systems Analysis Lab
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the % % "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, % distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to % the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.
%---------------------------------------------------------------------------------------------------------------------

%Create a simulation realization of a 3-dimensional (M=3) time-varying
%model of order p=2 and N=1000 samples.

function [Y,S,e]=createsimulationset(hetflag)

%TVcoef.mat contains a cell named coefcell. 
%coef is a 1xN cell. Each cell represents the MxM*p coefficient matrix 
%[A1(n) A2(n)] at each time point n (see Eq.23-24). 
%These matrices were created to produce a stable 3-dimensional (M=3) MVAR 
%model of order p=2 and N=1000 samples.

fprintf('\n--------- Creating Simulation... ')
tic;
load('TVcoef.mat')
M=3;
p=2;
N=1000;

if(hetflag==0)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Homoskedastic case
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %The variance of the noise that drives the MVAR model is 1 for all 3
    %time series (i.e. S=[1 0 0;0 1 0;0 0 1])
    Y=zeros(M,N);
    rng('shuffle') 
    e=randn(M,N);
    e=zscore(e')';
    S=[1 0 0;0 1 0;0 0 1];

    for k=p+1:N
        temp=Y(:,k-1:-1:k-p);
        Y(:,k)= coefcell{k}*temp(:)+e(:,k);
    end
elseif(hetflag==1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Heteroskedastic case
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %The variance of the noise that drives the MVAR model changes at time
    %point 551 and the noise covariance becomes S=[0.5 0 0;0 0.8 0;0 0 0.2]
    Y=zeros(M,N);
    rng('shuffle') 
    e=randn(M,N);
    e=zscore(e')';
    e(1,551:1000)=sqrt(0.5)*e(1,551:1000); 
    e(2,551:1000)=sqrt(0.8)*e(2,551:1000);
    e(3,551:1000)=sqrt(0.2)*e(3,551:1000);
    for k=1:550
        S{k}=[1 0 0;0 1 0;0 0 1];
    end
    for k=551:N
        S{k}=[0.5 0 0;0 0.8 0;0 0 0.2];
    end
    for k=p+1:N
        temp=Y(:,k-1:-1:k-p);
        Y(:,k)= coefcell{k}*temp(:)+e(:,k);
    end
end
fprintf('Total time: %f seconds\n',toc) 



