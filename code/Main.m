%---------------------------------------------------------------------------------------------------------------------
% Copyright (c) 2019 Kyriaki Kostoglou
% PhD Supervisor: Georgios Mitsis, Associate Professor, Bioengineering Department, McGill University, Montreal, Canada
% Biosignals and Systems Analysis Lab
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the % % "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, % distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to % the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.
%---------------------------------------------------------------------------------------------------------------------

close all;clear all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set parameter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param.pmax=10;           %Maximum MVAR model order to be considered when optimizing with the GA
param.pmaxGARCH=5;       %Maximum order to be considered when fitting GARCH models on the TV-MVAR residuals
param.metric=1;          %GA Fitness function (1 for multivariate AIC, 2 for multivariate BIC)
param.ignore=100;        %Number of time points to ignore due to initialization of the estimator (i.e. conventional/proposed KF)
param.smoothflag=1;      %Set to 1 to apply smoothing on the estimated TV-MVAR coefficients
param.hetflag=1;         %Set to 1 for the heteroskedastic case - The TV covariance of the MVAR residuals is estimated using GARCH models
param.measflag=1;        %Set to 1 to estimate the TV-MVAR measures (i.e. COH,PCOH,DC,gPDC) based on the obtained TV-MVAR coefficients
param.nfft=256;          %Number of points for the calculation of the TV-MVAR measures in the frequency domain 
param.fs=1;              %Sampling Frequency
siglab={'y1','y2','y3'}; %Labels for the time-series
param.M=3;               %Number of time-series (i.e. the dimension of the MVAR model)
heterosk=1;              %Set to 1 to simulate a TV-MVAR process driven by heteroskedastic noise or 0 for the homoskedastic case
                         %In the heteroskedastic case The variance of the driving noise is initially 1 (i.e. S=[1 0 0;0 1 0;0 0 1])                   
                         %for all time series but then the covariance changes to S=[0.5 0 0;0 0.8 0;0 0 0.2] at time point 551
                         %In the homoskedastic case the covariance is constant [1 0 0;0 1 0;0 0 1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create a simulation realization of a 3-dimensional (M=3) time-varying MVAR
%process of order p=2 and N=1000 samples. The model parameters were created
%in such a way that the model is stable at each time point. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Y,S,e]=createsimulationset(heterosk);
%Y:  M-dimensional TV-MVAR process of order p=2
%S:  True driving noise covariance
%e:  True driving noise
e=e(:,param.ignore:end);     %Initial samples are ignored so this vector can match the one returned from the MVAR analysis
                             %During recursive estimation ignore samples are ignored due to initialization 
                             
if(param.hetflag==1)         %If hetflag==1 then S is a cell.Each cell represents the covariance of the driving noise at each time point
   S=S(param.ignore:end);    %Again we ignore the first ignore cells
end

preal=2;                     %True TV-MVAR model order

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set GA options 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nworkers=4;            %Choose based on your computer specs
parpool(nworkers)      %Initialize parallel pool with nworkers
param.ga_opts = gaoptimset('TolFun',1e-12,'StallGenLimit',30,'Generations',100,'Display','iter','UseParallel','always');
% %To set parallel processing off use instead the following command:
%param.ga_opts = gaoptimset('TolFun',1e-12,'StallGenLimit',30,'Generations',100,'Display','iter');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GA optimization and estimation of the TV-MVAR model using the conventional KF 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[XKF,JKF]=GA_MVAR_KF(Y,param);
%Estimate the TV-MVAR model using the optimized GA parameters and save results to structure resultsKF
resultsKF=simulate_MVAR_KF(XKF,Y,param);                                                            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GA optimization and estimation of the TV-MVAR model using the proposed KF (MKFA) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[XMKFA, JMKFA]=GA_MVAR_MKFA(Y,param);
%Estimate the TV-MVAR model using the optimized GA parameters and save results
resultsMKFA=simulate_MVAR_MKFA(XMKFA,Y,param);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GA optimization and estimation of the TV-MVAR model using the proposed KF (MKFA) - improved version 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[XMKFA2, JMKFA2]=GA_MVAR_MKFA2(Y,param);
%Estimate the TV-MVAR model using the optimized GA parameters and save results
resultsMKFA2=simulate_MVAR_MKFA2(XMKFA2,Y,param);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the real and estimated TV-MVAR measures
if(param.measflag==1)
    load('TVcoef.mat')
    meth={'real','KF','MKFA','MKFA2'};
    %Estimate the real TV-MVAR measures
    real_MVARmeasures=estimate_MVARmeasures(coefcell(param.ignore:end),S,param.nfft,param.fs,param.M,preal,param.hetflag);
    %Plot the real TV-MVAR measures
    plot_MVARmeasures(real_MVARmeasures,siglab,meth{1})
    %Plot the estimated TV-MVAR measures based on the conventional KF
    plot_MVARmeasures(resultsKF.MVARmeasures,siglab,meth{2})
    %Plot the estimated TV-MVAR measures based on the proposed KF (MKFA)    
    plot_MVARmeasures(resultsMKFA.MVARmeasures,siglab,meth{3})
    %Plot the estimated TV-MVAR measures based on the proposed KF (MKFA) - improved version   
    plot_MVARmeasures(resultsMKFA2.MVARmeasures,siglab,meth{4})    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimate NMSE between real and estimated TV-MVAR measures (see Eq.26)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for mT=1:param.M
    for mD=1:param.M
        %Conventional KF
        eCOHKF(mT,mD)=compute_NMSE(resultsKF.MVARmeasures.COH{mT,mD},real_MVARmeasures.COH{mT,mD});
        ePCOHKF(mT,mD)=compute_NMSE(resultsKF.MVARmeasures.PCOH{mT,mD},real_MVARmeasures.PCOH{mT,mD});   
        eDCKF(mT,mD)=compute_NMSE(resultsKF.MVARmeasures.DC{mT,mD},real_MVARmeasures.DC{mT,mD});   
        egPDCKF(mT,mD)=compute_NMSE(resultsKF.MVARmeasures.gPDC{mT,mD},real_MVARmeasures.gPDC{mT,mD});   

        %Proposed KF (MKFA)
        eCOHMKFA(mT,mD)=compute_NMSE(resultsMKFA.MVARmeasures.COH{mT,mD},real_MVARmeasures.COH{mT,mD});
        ePCOHMKFA(mT,mD)=compute_NMSE(resultsMKFA.MVARmeasures.PCOH{mT,mD},real_MVARmeasures.PCOH{mT,mD});   
        eDCMKFA(mT,mD)=compute_NMSE(resultsMKFA.MVARmeasures.DC{mT,mD},real_MVARmeasures.DC{mT,mD});   
        egPDCMKFA(mT,mD)=compute_NMSE(resultsMKFA.MVARmeasures.gPDC{mT,mD},real_MVARmeasures.gPDC{mT,mD});   
        
        %Proposed KF (MKFA2) - improved version
        eCOHMKFA2(mT,mD)=compute_NMSE(resultsMKFA2.MVARmeasures.COH{mT,mD},real_MVARmeasures.COH{mT,mD});
        ePCOHMKFA2(mT,mD)=compute_NMSE(resultsMKFA2.MVARmeasures.PCOH{mT,mD},real_MVARmeasures.PCOH{mT,mD});   
        eDCMKFA2(mT,mD)=compute_NMSE(resultsMKFA2.MVARmeasures.DC{mT,mD},real_MVARmeasures.DC{mT,mD});   
        egPDCMKFA2(mT,mD)=compute_NMSE(resultsMKFA2.MVARmeasures.gPDC{mT,mD},real_MVARmeasures.gPDC{mT,mD});        
    end
end

%NMSE for TV-MVAR measures using the conventional KF
NMSE_COHKF=mean(eCOHKF(:))*100;
NMSE_PCOHKF=mean(ePCOHKF(:))*100;
NMSE_DCKF=mean(eDCKF(:))*100;
NMSE_PDCKF=mean(egPDCKF(:))*100;

%NMSE for TV-MVAR measures using the proposed KF (MKFA)
NMSE_COHMKFA=mean(eCOHMKFA(:))*100;
NMSE_PCOHMKFA=mean(ePCOHMKFA(:))*100;
NMSE_DCMKFA=mean(eDCMKFA(:))*100;
NMSE_PDCMKFA=mean(egPDCMKFA(:))*100;

%NMSE for TV-MVAR measures using the proposed KF (MKFA2) - improved version
NMSE_COHMKFA2=mean(eCOHMKFA2(:))*100;
NMSE_PCOHMKFA2=mean(ePCOHMKFA2(:))*100;
NMSE_DCMKFA2=mean(eDCMKFA2(:))*100;
NMSE_PDCMKFA2=mean(egPDCMKFA2(:))*100;

%Create table and print on command window
fprintf('\n----------------------')
fprintf('\nTV-MVAR measures NMSE ')
fprintf('\n----------------------')
Measure={'COH';'PCOH';'DC';'gPDC'};
NMSE_KF=[NMSE_COHKF;NMSE_PCOHKF;NMSE_DCKF;NMSE_PDCKF];
NMSE_MKFA=[NMSE_COHMKFA;NMSE_PCOHMKFA;NMSE_DCMKFA;NMSE_PDCMKFA];
NMSE_MKFA2=[NMSE_COHMKFA2;NMSE_PCOHMKFA2;NMSE_DCMKFA2;NMSE_PDCMKFA2];
T = table(Measure,NMSE_KF,NMSE_MKFA,NMSE_MKFA2,'VariableNames',{'Measure','KF','MKFA','MKFA2'})



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimate MSE between real and estimated TV-MVAR coefficients (see Eq.26)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(param.smoothflag==1) %if smoothflag=1 then we use the smoothed vector of estimated coefficients
    %Conventional KF
    MSEcKF=compute_NMSE(resultsKF.thvecsmooth,coefvec(:,param.ignore:end));   

    %Proposed KF (MKFA)
    MSEcMKFA=compute_NMSE(resultsMKFA.thvecsmooth,coefvec(:,param.ignore:end));    

    %Proposed KF (MKFA2) - improved version
    MSEcMKFA2=compute_NMSE(resultsMKFA2.thvecsmooth,coefvec(:,param.ignore:end));    
else
    %Conventional KF
    MSEcKF=compute_NMSE(resultsKF.thvec,coefvec(:,param.ignore:end));   

    %Proposed KF (MKFA)
    MSEcMKFA=compute_NMSE(resultsMKFA.thvec,coefvec(:,param.ignore:end));    

    %Proposed KF (MKFA2) - improved version
    MSEcMKFA2=compute_NMSE(resultsMKFA2.thvec,coefvec(:,param.ignore:end));   
end


%Create table and print on command window
fprintf('\n----------------------------------------------------')
fprintf('\nMSE between true and predicted TV-MVAR coefficients ')
fprintf('\n----------------------------------------------------')
method={'KF';'MKFA';'MKFA2'};
MSEc=[MSEcKF;MSEcMKFA;MSEcMKFA2];
T = table(method,MSEc)



save('RES.mat')




