%---------------------------------------------------------------------------------------------------------------------
% Copyright (c) 2019 Kyriaki Kostoglou
% PhD Supervisor: Georgios Mitsis, Associate Professor, Bioengineering Department, McGill University, Montreal, Canada
% Biosignals and Systems Analysis Lab
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the % % "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, % distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to % the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.
%---------------------------------------------------------------------------------------------------------------------

%Estimate the TV-MVAR model using the optimized GA parameters and the conventional KF as estimator
%If smoothflag is 1, then smoothing is applied on the estimated TV-MVAR coefficients
%If hetflag is 1, then it is assumed that the TV-MVAR residuals are heteroskedastic and GARCH models are used to fit the variance of the %residuals
%If measflag is 1, then the TV-MVAR measures of COH, PCOH, DC, gPDC are estimated in the frequency domain

function results=simulate_MVAR_KF(X,Y,param)

metric=param.metric;
ignore=param.ignore;
pmaxGARCH=param.pmaxGARCH;
smoothflag=param.smoothflag;
hetflag=param.hetflag;
measflag=param.measflag;
fs=param.fs;
nfft=param.nfft;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X:          row vector that contains the hyperparameters optimized by the GA
% Y:          MxN matrix (each time series is placed in a row)
% metric:     fitness function (1 for multivariate AIC, 2 for multivariate BIC)
% ignore:     number of samples to ignore due to initialization
% pmaxGARCH:  maximum order to be considered when fitting GARCH models on the TV-MVAR residuals
% smoothflag: set to 1 for applying smoothing to the estimated parameters
% hetflag:    set to 1 to fit GARCH models on the TV-MVAR residuals
% measflag:   set to 1 to estimate the TV-MVAR measures (i.e. COH,PCOH,DC,gPDC) based on the obtained TV-MVAR coefficients
% fs:         sampling frequency
% nfft:       number of points for the calculation of the TV-MVAR measures in the frequency domain 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p=X(1);                 %Model order for the TV-MVAR model
[M,N]=size(Y);    
Mp=M*p;    
totPar=M*Mp;            %Total number of model parameters M*M*p
pin=X(end);             %P0: Initial value for the covariance matrix P as optimized by the GA
P = pin*eye(Mp);        %P:  Parameter estimation error covariance matrix
R2=X(2);                %R2: Measurement noise variance as optimized by the GA
R1=X(3)*eye(Mp);        %R1: Diagonal process noise covariance matrix as optimized by the GA
yhat=zeros(M,N);        %Predicted output
th=zeros(M,Mp);         %MxMp Initial parameter vector (zero)     
e=zeros(M,N);
thvec=zeros(totPar,N);  %Matrix were the TV coefficients are saved 

%Save hyperparameters to structure results
results.X=X;            %GA optimum solution
results.p=p;            %GA optimum model order
results.pin=pin;        %GA optimum P0
results.R2=R2;          %GA optimum R2
results.R1=R1;          %GA optimum R1
results.M=M;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Conventional Kalman Filter (see Eqs.1-7)
fprintf('\n->Estimating the time-varying MVAR model... ')  
tic;
for k=p+1:N
    temp= Y(:,k-1:-1:k-p);
    phi=temp(:);
    phit=phi';
    yh=th*phi;
    e(:,k)=Y(:,k)-yh;
    Pprior{k}=P+R1;       %needed for smoothing the estimated coefficients
    phitP=phit*P;
    rt=phitP*phi;
    K=(P*phi)/(R2+rt); 
    P=P+R1-K*phitP;
    Ppost{k}=P;           %needed for smoothing the estimated coefficients
    th=th+(K*e(:,k)')';  
    yhat(:,k)=yh;
    thcell{k}=th;         %the estimated coefficients are saved in cells and matrices (for ease of use in smoothing and visualization)
    thvec(:,k)=th(:);
end
fprintf('Total time: %f',toc)  


%In structure results all signals are saved from time point ignore and after. ignore samples are ignored due to estimator initialization
results.thcell=thcell(ignore:end);     %estimated TV-MVAR coefficients in cell form (each cell is on time point and size MxM*p)
results.thvec=thvec(:,ignore:end);     %estimated TV-MVAR coefficients in matrix form (each column is one time point-all coefficients are   
                                       %concantenated in a column of length M*M*px1)  
results.YHAT=yhat(:,ignore:end);       %predicted Y
results.YREAL=Y(:,ignore:end);         %Y
results.e=e(:,ignore:end);             %error

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
        J=NN*log(ds)+2*totPar;        %Multivariate AIC
    elseif(metric==2)
        J=NN*log(ds)+log(NN)*totPar;  %Multivariate BIC
    end
end
results.J=J;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Smoothing of the TV estimated parameters (see Eq.13-14)
if(smoothflag==1)
   fprintf('\n->Applying smoothing to the estimated MVAR parameters... ')  
   tic;
   [thcellsmooth, thvecsmooth, yhsmooth]=smooth_MVAR_KF(Pprior,Ppost,Y,thcell,M,N,p,ignore);
   %thcellsmooth: smoothed estimated TV-MVAR coefficients saved in cells (each cell is one time point)
   %thvecsmooth:  smoothed estimated TV-MVAR coefficients saved in matrices (each column is one time point)
   %yhsmooth:     predicted MxN-ignore output based on the estimated smoothed TV-MVAR coefficients
   results.thcellsmooth=thcellsmooth;
   results.thvecsmooth=thvecsmooth;
   results.YHATsmooth=yhsmooth;
   fprintf('Total time: %f',toc)  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GARCH Analysis if hetflag is 1
if(hetflag==1)
   fprintf('\n->Estimating the time-varying GARCH covariance for the MVAR residuals...\n')   
   tic;
   if(smoothflag==1)
      e=Y-yhsmooth;
   end
   resGARCH=estimate_GARCH(e,metric,ignore,pmaxGARCH,S);%the residuals are assumed to be heteroskedastic and the GARCH approach is applied
   results.GARCH=resGARCH;                              %save the results of the GARCH analysis to structure results
   S=resGARCH.SGARCH;                                   %GARCH estimated residual TV covariance (in cells - each cell is one time point)
   fprintf('\nTotal time: %f',toc)   
end

results.S=S;    %the covariance of the residuals (if hetflag=1 then this is TV and each cell is one time point. 
                %Each cell is a matrix of size MxM representing the covariance of the residuals at one time point) 
                %If hetflag is not 1 then S is simply cov(e(:,ignore:end))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimate TV-MVAR measures (Eqs.15-20) if measflag==1
if(measflag==1)
   fprintf('\n\n->Estimating the time-varying MVAR measures... ')   
   tic;
   if(smoothflag==1) 
        %If smoothflag=1 use the smoothed estimated TV-MVAR coefficients
        results.MVARmeasures=estimate_MVARmeasures(thcellsmooth,S,nfft,fs,M,p,hetflag);
   else
        %else use the unsmoothed estimated TV-MVAR coefficients
        results.MVARmeasures=estimate_MVARmeasures(thcell,S,nfft,fs,M,p,hetflag);
   end
   fprintf('Total time: %f',toc)   
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Structure results contains the following:
%results.X=X;                          %GA optimum solution
%results.p=p;                          %GA optimum model order
%results.pin=pin;                      %GA optimum P0
%results.R2=R2;                        %GA optimum R2
%results.R1=R1;                        %GA optimum R1
%results.M=M;                          %Number of MVAR time-series
%results.thcell=thcell(ignore:end);    %estimated TV-MVAR coefficients in cell form (each cell is on time point and size MxM*p)
%results.thvec=thvec(:,ignore:end);    %estimated TV-MVAR coefficients in matrix form (each column is one time point-all coefficients are   
                                       %concantenated in a column of length M*M*px1)    
                                        
%results.YHAT=yhat(:,ignore:end);      %predicted Y
%results.YREAL=Y(:,ignore:end);        %Y
%results.e=e(:,ignore:end);            %error
%results.J=J;                          %minimum fitness function

%-if smoothflag is 1 then also the following are included:
%results.thcellsmooth=thcellsmooth;    %smoothed estimated TV-MVAR coefficients in cell form
%results.thvecsmooth=thvecsmooth;      %smoothed estimated TV-MVAR coefficients in matrix form
%results.YHATsmooth=yhsmooth;          %predicted Y using the smoothed estimated TV-MVAR coefficients

%-if hetflag is 1 then also the following are included:
%results.GARCH=resGARCH;               %save the results of the GARCH analysis to structure results 

%-if measflag is 1 then compute TV-MVAR measures of COH, COH, DC, gPDC
%results.MVARmeasures

%-also:
%results.S=S;                          %the covariance of the residuals (if hetflag=1 then this is TV and each cell is one time point. 
                                       %Each cell is a matrix of size MxM representing the covariance of the residuals at one time point) 
                                       %If hetflag is not 1 then S is simply cov(e(:,ignore:end))