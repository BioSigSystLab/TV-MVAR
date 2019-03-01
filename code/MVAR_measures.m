%---------------------------------------------------------------------------------------------------------------------
% Copyright (c) 2019 Kyriaki Kostoglou
% PhD Supervisor: Georgios Mitsis, Associate Professor, Bioengineering Department, McGill University, Montreal, Canada
% Biosignals and Systems Analysis Lab
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the % % "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, % distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to % the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.
%---------------------------------------------------------------------------------------------------------------------

%Estimate MVAR measures (i.e. COH,PCOH,DC,gPDC) 

function results = MVAR_measures(th,Se,nfft,fs,M,p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% th:    Estimated MVAR coefficients in a matrix form (MxM*p)
% Se:    Estimated MVAR residual covariance matrix (MxM)
% nfft:  Number of points for the calculation of the MVAR measures (i.e. COH,PCOH,DC,gPDC) in the frequency domain 
% fs:    Sampling Frequency
% M:     Number of time-series
% p:     MVAR model order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
for k=1:p
    A{k}=th(:,(k-1)*M+1:(k-1)*M+M);             %Extracting the k-th MxM block from the estimated MVAR coefficient matrix
end
f = (0:nfft-1)*(fs/(2*nfft));                   %Frequency axis
dSe1=repmat(diag(sqrt(Se))',M,1);               %Needed for DC estimation                          
idSe1=1./repmat(diag(sqrt(Se)),1,M);            %Needed for gPDC estimation 
dSe2=repmat(diag(Se)',M,1);                     %Needed for DC estimation   
idSe2=1./repmat(diag(Se),1,M);                  %Needed for gPDC estimation 
h=cell(nfft,1);
sp=cell(nfft,1);
pp=cell(nfft,1);
coh=cell(nfft,1);
pcoh=cell(nfft,1);
dc=cell(nfft,1);
gpdc=cell(nfft,1);

for n=1:nfft 
    Af=zeros(M,M);
    for k=1:p
        Af=Af+A{k}*exp(-1i*2*pi*f(n)*k/fs);                     %Coefficient Matrix in frequency domain (see Section II.E 1st paragraph)
    end
    Aft=eye(M)-Af;
    h{n}=inv(Aft);                                                        %Tranfer Function Matrix
    sp{n}=h{n}*Se*h{n}';                                                  %Spectral Matrix (Eq.15)
    pp{n}=inv(sp{n});                                                     %Inverse Spectral Matrix (Eq.16)
    coh{n}=sp{n}./sqrt(diag(sp{n})*diag(sp{n})');                         %Coherence (Eq.17)    
    pcoh{n}=pp{n}./sqrt(diag(pp{n})*diag(pp{n})');                        %Partial Coherence (Eq.19)
    dc{n}=(dSe1.*h{n})./sqrt(repmat(sum(dSe2.*(abs(h{n}).^2),2),1,M));    %Directed Coherence (Eq.18)
    gpdc{n}=(idSe1.*Aft)./sqrt(repmat(sum(idSe2.*(abs(Aft).^2),1),M,1));  %Generalized Partial Directed Coherence (Eq.20)) 
end

H=cell(M);
S=cell(M);
P=cell(M);
COH=cell(M);
PCOH=cell(M);
DC=cell(M);
gPDC=cell(M);


% For ease of use the estimated TV-MVAR measures are saved on cells with
% indexes {mT,mD}, where mT is the target time-series and mD is the driver
% time-series.

for mT=1:M
    for mD=1:M
        H{mT,mD}=abs(cellfun(@(X)X(mT,mD),h));                     %taking the modulus of H
        S{mT,mD}=abs(cellfun(@(X)X(mT,mD),sp));                    %taking the modulus of S
        S{mT,mD}=abs(S{mT,mD}/max(S{mT,mD}));                      %normalization of S
        P{mT,mD}=abs(cellfun(@(X)X(mT,mD),pp));                    %taking the modulus of P
        COH{mT,mD}=abs(cellfun(@(X)X(mT,mD),coh)).^2;              %taking the squared modulus of COH
        PCOH{mT,mD}=abs(cellfun(@(X)X(mT,mD),pcoh)).^2;            %taking the squared modulus of PCOH
        DC{mT,mD}=abs(cellfun(@(X)X(mT,mD),dc)).^2;                %taking the squared modulus of DC
        gPDC{mT,mD}=abs(cellfun(@(X)X(mT,mD),gpdc)).^2;            %taking the squared modulus of gPDC
    end
end

%Save estimated MVAR measures to structure results
results.H=H;            
results.S=S;
results.P=P;
results.COH=COH;
results.PCOH=PCOH;
results.DC=DC;
results.gPDC=gPDC;
results.f=f;
results.fs=fs;
results.nfft=nfft;