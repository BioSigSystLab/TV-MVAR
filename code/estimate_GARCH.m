%---------------------------------------------------------------------------------------------------------------------
% Copyright (c) 2019 Kyriaki Kostoglou
% PhD Supervisor: Georgios Mitsis, Associate Professor, Bioengineering Department, McGill University, Montreal, Canada
% Biosignals and Systems Analysis Lab
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the % % "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, % distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to % the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.
%---------------------------------------------------------------------------------------------------------------------

%Estimate TV covariance of the MVAR residuals using GARCH models (see Eq.22)

function results=estimate_GARCH(e,metric,ignore,pmax,S)
M=size(e,1);
for m=1:M  
    %Testing different GARCH model orders for time-series m
    count=0;
    md=[];
    aic=[];
    bic=[];
    for q=0:pmax
        for r=1:pmax
            flag=0;
            fprintf('\nFitting GARCH to time-series %d - Testing model orders: [q=%d,r=%d] ',m,q,r)
            Mdl = garch(q,r);
            try
                %The following command may return the following
                %error:"Non-zero degree P requires a non-zero degree Q." is
                %there is no underlying true heteroskedasticity. Therefore
                %we try to catch the error.
                [EstMdl,EstParamCov,logL]  = estimate(Mdl,e(m,:)','Display','off');
                %If it finally goes through then flag is set to 1
                flag=1;
            catch
                %do nothing
            end
            if(flag==1)
                count=count+1;
                md(count,:)=[q r];
                [aic(count),bic(count)] = aicbic(logL,q+r,length(e(m,:)));   %compute AIC/BIC scores based on the estimated GARCH model         
            end
        end
    end
    
    %if aic or bic are nonempty vectors then we select the optimal model order based on the aic or bic scores
    %if aic or bic are empty vectors then it means that there is no real
    %heteroskedasticity in the data to be GARCH modeled 
    if(isempty(aic)==0)       
        if(metric==1)
            [mm,ii]=min(aic);
        else
            [mm,ii]=min(bic);
        end
        %Find optimal GARCH model orders by selecting the pair [q,r] with the minimum AIC/BIC score 
        fprintf('\nOptimal model order found: [q=%d,r=%d]\n',md(ii,:))
        Mdl = garch(md(ii,1),md(ii,2));
        [EstMdl]  = estimate(Mdl,e(m,:)');
        temp=infer(EstMdl,e(m,:)');
        results.EstMdl{m}=EstMdl;
        %Save TV variance as predicted by the fitted GARCH model for each residual time-series
        results.yGARCH(m,:)=smooth(temp(ignore:end),100);   %smooth the garch variance prediction
        results.hete(m)=1;
    else
        fprintf('\nNo heteroskedasticity detected')
        results.hete(m)=0;
    end
        
end

%Create diagonal TV covariance of the MVAR residuals
newN=size(results.yGARCH,2);
for k=1:newN
    results.SGARCH{k}=zeros(M,M);
    for m=1:M
        if(results.hete(m)==0)
            %if GARCH modeling failed, it means that there is no
            %heteroskedasticity in the residuals and we just use the
            %estimated variance on the whole residual time-series
            results.SGARCH{k}(m,m)=S(m,m);
        else
            %else we use the predicted GARCH TV variance
            results.SGARCH{k}(m,m)=results.yGARCH(m,k);
        end
    end         
end
end



%results is a structure that contains the following:
%results.EstMdl{m}=EstMdl;         %optimum GARCH estimated model for the mth time-series residuals
%results.yGARCH(m,:)               %predicted TV GARCH variance for the mth time-series residuals
%results.SGARCH{k}                 %TV GARCH covariance of the residuals at time point k
%results.hete(m)                   %0 if there is no underlying heteroskedasticity in the data else it is set to 1 and heteroskedasticity was modeled using GARCH models 