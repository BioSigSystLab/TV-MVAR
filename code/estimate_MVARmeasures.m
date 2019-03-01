%---------------------------------------------------------------------------------------------------------------------
% Copyright (c) 2019 Kyriaki Kostoglou
% PhD Supervisor: Georgios Mitsis, Associate Professor, Bioengineering Department, McGill University, Montreal, Canada
% Biosignals and Systems Analysis Lab
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the % % "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, % distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to % the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.
%---------------------------------------------------------------------------------------------------------------------

%Estimate TV-MVAR measures (i.e. COH,PCOH,DC,gPDC) 

function results=estimate_MVARmeasures(th,S,nfft,fs,M,p,hetflag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% th:      Estimated TV-MVAR coefficients in a cell of length N. Each cell contains the MxM*p matrix of the estimated MVAR coefficients at %          each time point
% S:       Estimated MVAR residual covariance matrix. If hetflag is 1 then S is a cell of length N. Each cell contains the GARCH estimated  
%          covariance at each time point. Otherwise it is a MxM matrix.
% nfft:    Number of points for the calculation of the MVAR measures (i.e. COH,PCOH,DC,gPDC) in the frequency domain 
% fs:      Sampling Frequency
% M:       Number of time-series
% p:       MVAR model order
% hetflag: Set to 1 for the heteroskedastic case (i.e. the covariance of the MVAR residuals is estimated using GARCH models)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Estimate MVAR measures at each time point
for k=1:length(th)
    if(hetflag==1)
       %if hetflag==1 then S is TV (cell of length N - as estimated by the GARCH models) and therefore at each time point we need to compute the MVAR measures based on the obtained TV MVAR coefficients and the estimated GARCH covariance (S) of the model residuals 
       res{k}=MVAR_measures(th{k},S{k},nfft,fs,M,p);
    else
       %else S is just an MxM matrix
       res{k}=MVAR_measures(th{k},S,nfft,fs,M,p);
    end
end

%The TV-MVAR measures are saved on cells with indexes {mT,mD}, where mT is the target time-series and mD is the driver time-series.
%Each column of a cell represents time. E.g. results.PCOH{1,2}(:,10) is the PCOH from time-series 2 to time-series 1 at time point 10. 
for k=1:length(th)
     for mT=1:M        
        for mD=1:M
            results.SP{mT,mD}(:,k)=res{k}.S{mT,mD};
            results.P{mT,mD}(:,k)=res{k}.P{mT,mD};
            results.COH{mT,mD}(:,k)=res{k}.COH{mT,mD};
            results.PCOH{mT,mD}(:,k)=res{k}.PCOH{mT,mD};            
            results.DC{mT,mD}(:,k)=res{k}.DC{mT,mD};
            results.gPDC{mT,mD}(:,k)=res{k}.gPDC{mT,mD};
            results.f=res{k}.f;
        end
     end
end
results.fs=fs;
results.M=M;



%results is a structure that contains the following:
%results.SP{mT,mD}(:,k)      %spectral power density matrix from time-series mD to time-series mT at time point k
%results.P{mT,mD}(:,k)       %inverse spectral power density matrix from time-series mD to time-series mT at time point k
%results.COH{mT,mD}(:,k)     %Coherence from time-series mD to time-series mT at time point k
%results.PCOH{mT,mD}(:,k)    %Parital Coherence from time-series mD to time-series mT at time point k      
%results.DC{mT,mD}(:,k)      %Directed Coherence from time-series mD to time-series mT at time point k
%results.gPDC{mT,mD}(:,k)    %Generalized Partial Directed Coherence from time-series mD to time-series mT at time point k
%results.f                   %frequency axis based on nfft and fs
%results.fs                  %sampling frequency
%results.M                   %number of MVAR time-series
%