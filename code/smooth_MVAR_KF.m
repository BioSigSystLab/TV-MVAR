%---------------------------------------------------------------------------------------------------------------------
% Copyright (c) 2019 Kyriaki Kostoglou
% PhD Supervisor: Georgios Mitsis, Associate Professor, Bioengineering Department, McGill University, Montreal, Canada
% Biosignals and Systems Analysis Lab
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the % % "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, % distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to % the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.
%---------------------------------------------------------------------------------------------------------------------

%Smoothing of the estimated TV MVAR_KF parameters using Rauch-Tung-Striebel
%fixed-interval equations (see Eq.13-14)

function [thcellsmooth, thvecsmooth, yhsmooth]=smooth_MVAR_KF(Pprior,Ppost,Y,thcell,M,N,p,ignore)

thvecsmooth=zeros(M*M*p,N);
yhsmooth=zeros(M,N);
for k=N:-1:p+1
    if(k==N)
        thcellsmooth{N}=thcell{N};
        tempc=thcell{k};
        thvecsmooth(:,k)=tempc(:);
    else
        W=Ppost{k}/Pprior{k+1};
        thcellsmooth{k}=thcell{k}+(thcellsmooth{k+1}-thcell{k})*W;
        tempc=thcellsmooth{k};
        thvecsmooth(:,k)=tempc(:);      
    end
    temp= Y(:,k-1:-1:k-p);
    phi=temp(:);
    yhsmooth(:,k)=thcellsmooth{k}*phi;      
end
thcellsmooth=thcellsmooth(ignore:end);
thvecsmooth=thvecsmooth(:,ignore:end);
end


