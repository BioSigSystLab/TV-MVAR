%---------------------------------------------------------------------------------------------------------------------
% Copyright (c) 2019 Kyriaki Kostoglou
% PhD Supervisor: Georgios Mitsis, Associate Professor, Bioengineering Department, McGill University, Montreal, Canada
% Biosignals and Systems Analysis Lab
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the % % "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, % distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to % the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.
%---------------------------------------------------------------------------------------------------------------------

%Smoothing of the estimated TV MVAR_MKFA parameters using Rauch-Tung-Striebel
%fixed-interval equations (see Eq.13-14)

function [thcellsmooth, thvecsmooth, yhsmooth]=smooth_MVAR_MKFA(Pprior,Ppost,Y,thvec2,M,N,p,ignore)

thvecsmooth2=zeros(M*M*p,N);
thvecsmooth=zeros(M*M*p,N);
yhsmooth=zeros(M,N);

for k=N:-1:p+1
    if(k==N)
        thvecsmooth2(:,k)=thvec2(:,N);
        temp=reshape(thvecsmooth2(:,k),M*p,M)';
        thvecsmooth(:,k)=temp(:);
        thcellsmooth{k}=temp;
    else
        W=Ppost{k}/Pprior{k+1};
        thvecsmooth2(:,k)=thvec2(:,k)+W*(thvecsmooth2(:,k+1)-thvec2(:,k));
        temp=reshape(thvecsmooth2(:,k),M*p,M)';
        thvecsmooth(:,k)=temp(:);
        thcellsmooth{k}=temp;
    end
    temp=Y(:,k-1:-1:k-p);
    temp=temp(:);
    for m=1:M
        phit(m,(m-1)*M*p+1:(m-1)*M*p+M*p)=temp;
    end
    yhsmooth(:,k)=phit*thvecsmooth2(:,k);    
end
thcellsmooth=thcellsmooth(ignore:end);
thvecsmooth=thvecsmooth(:,ignore:end);
end


