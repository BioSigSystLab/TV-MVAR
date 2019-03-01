%---------------------------------------------------------------------------------------------------------------------
% Copyright (c) 2019 Kyriaki Kostoglou
% PhD Supervisor: Georgios Mitsis, Associate Professor, Bioengineering Department, McGill University, Montreal, Canada
% Biosignals and Systems Analysis Lab
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the % % "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, % distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to % the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.
%---------------------------------------------------------------------------------------------------------------------

%Function that initiates the GA optimization of the MVAR model hyperparameters based on the conventional KF 

function [XKF,JKF]=GA_MVAR_KF(Y,param)

metric=param.metric;
ignore=param.ignore;
ga_opts=param.ga_opts;
pmax=param.pmax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y:          MxN matrix (each time series is placed in a row)
% metric:     Fitness function (1 for multivariate AIC, 2 for multivariate BIC)
% ignore:     Number of samples to ignore due to initialization of the estimators (i.e. conventional/proposed KF)
% ga_opts:    GA options
% pmax:       Maximum MVAR model order to be considered when optimizing with the GA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n-------- Initiating GA MVAR model optimization based on the conventional KF --------\n\n') 
h = @(X) MVAR_KF(X,Y,metric,ignore);                                         
%In case you want to use mex files for faster runtime type instead: h = @(X) MVAR_KF_mex(X,Y,metric,ignore);
nvars=4;                                                        %Number of hyperparameters optimized by the GA (p,R1,R2,P0)
LB=[1 0 0 0];                                                   %Lower bound for the hyperparameters values 
UB=[pmax Inf Inf Inf];                                          %Upper bound for the hyperparameters values 
warning('off','all')
tic;
[XKF, JKF] = ga(h, nvars,[],[],[],[],LB,UB,[],[1],ga_opts);    %XKF: Optimized hyperparameters returned by the GA  
                                                               %JKF: Obtained minimum value of the fitness function                     
timeKF=toc;
fprintf('\n-------- Optimization Completed - Total time: %f seconds --------\n\n',timeKF) 

end