%-----------------------------------------------------------------------------------------------------------
% Copyright (c) 2019 Kyriaki Kostoglou
% PhD Supervisor: Georgios Mitsis, Associate Professor, Bioengineering Department, McGill University, Montreal, Canada
% Biosignals and Systems Analysis Lab
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the % % "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, % distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to % the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.
%-----------------------------------------------------------------------------------------------------------

%Compute Mean Squared Error (NMSE) between true and predicted TV-MVAR measures coefficients
%Note that the predicted model order must be the same as the true model order. Only then the total number of predicted coefficients M*M*p %is equal to the number of true coefficients (M*M*ptrue). In all realization the GA is able to detect the correct model order and %therefore this is not a problem most of the times.

function MSE=compute_MSE(MAT,MATr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MAT:   Predicted TV-MVAR coefficients; a matrix of size M*M*pxN-ignore+1 (each row is one coefficient)
%MATr:  True TV-MVAR coefficients; a matrix of size M*M*pxN-ignore+1 (each row is one coefficient)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MSE=mean(sum((MAT-MATr).^2,2)./size(MATr,2));