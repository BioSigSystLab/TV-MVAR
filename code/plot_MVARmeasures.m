%---------------------------------------------------------------------------------------------------------------------
% Copyright (c) 2019 Kyriaki Kostoglou
% PhD Supervisor: Georgios Mitsis, Associate Professor, Bioengineering Department, McGill University, Montreal, Canada
% Biosignals and Systems Analysis Lab
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the % % "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, % distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to % the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% The Software is provided "as is", without warranty of any kind.
%---------------------------------------------------------------------------------------------------------------------

%Plot TV-MVAR measures

function plot_MVARmeasures(res,siglab,meth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%res:     Structure returned by the function estimate_MVARmeasures 
%siglab:  Labels for the time-series 
%meth:    A string that indicates that the plotted TV-MVAR measures were estimated using a specific algorithm (i.e. conventional KF, proposed KF, true etc.) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f=res.f;        %Frequency axis
fs=res.fs;      %Sampling frequency
M=res.M;        %Number of time-series

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot COH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
[map, descriptorname, description] = colorcet('L3');
colormap(map)
t=(0:size(res.COH{1,1},2)-1)/fs;
for mT=1:M
    for mD=1:M
        if(mT==mD)
            subplot(M,M,(mT-1)*M+mD); 
            imagesc(t,f,res.SP{mT,mD})
            set(gca,'YDir','normal')
            title(sprintf('%s_S %s',meth,siglab{mT}));
            caxis([0 1])
        else
            subplot(M,M,(mT-1)*M+mD); 
            imagesc(t,f,res.COH{mT,mD})
            set(gca,'YDir','normal')
            title(sprintf('%s_C_O_H %s -> %s',meth,siglab{mD},siglab{mT}));
            caxis([0 1])
        end
        if(mD==1)
            ylabel('Hz')
        else
           set(gca, 'YTickLabelMode', 'Manual')
           set(gca, 'YTick', [])
        end
        if(mT==M)
            xlabel('n(sec)')
        else
           set(gca, 'XTickLabelMode', 'Manual')
           set(gca, 'XTick', [])
        end
        set(gca,'fontsize',12,'fontweight','bold')
    end
end
saveas(gcf,sprintf('%sCOH.png',meth))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot PCOH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
[map, descriptorname, description] = colorcet('L3');
colormap(map)
for mT=1:M
    for mD=1:M
        if(mT==mD)
            subplot(M,M,(mT-1)*M+mD); 
            imagesc(t,f,res.P{mT,mD})
            set(gca,'YDir','normal')
            title(sprintf('%s_P %s',meth,siglab{mT}));
            caxis([0 1])
        else
            subplot(M,M,(mT-1)*M+mD); 
            imagesc(t,f,res.PCOH{mT,mD})
            set(gca,'YDir','normal')
            title(sprintf('%s_P_C_O_H %s -> %s',meth,siglab{mD},siglab{mT}));
            caxis([0 1])
        end
        if(mD==1)
            ylabel('Hz')
        else
           set(gca, 'YTickLabelMode', 'Manual')
           set(gca, 'YTick', [])
        end
        if(mT==M)
            xlabel('n(sec)')
        else
           set(gca, 'XTickLabelMode', 'Manual')
           set(gca, 'XTick', [])
        end
        set(gca,'fontsize',12,'fontweight','bold')
    end
end
saveas(gcf,sprintf('%sPCOH.png',meth))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot DC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
[map, descriptorname, description] = colorcet('L3');
colormap(map)
for mT=1:M
    for mD=1:M
        subplot(M,M,(mT-1)*M+mD); 
        imagesc(t,f,res.DC{mT,mD})
        set(gca,'YDir','normal')
        title(sprintf('%s_D_C %s -> %s',meth,siglab{mD},siglab{mT}));
        caxis([0 1])
        if(mD==1)
            ylabel('Hz')
        else
           set(gca, 'YTickLabelMode', 'Manual')
           set(gca, 'YTick', [])
        end
        if(mT==M)
            xlabel('n(sec)')
        else
           set(gca, 'XTickLabelMode', 'Manual')
           set(gca, 'XTick', [])
        end
        set(gca,'fontsize',12,'fontweight','bold')
    end
end
saveas(gcf,sprintf('%sDC.png',meth))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot gPDC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
[map, descriptorname, description] = colorcet('L3');
colormap(map)
for mT=1:M
    for mD=1:M
        subplot(M,M,(mT-1)*M+mD); 
        imagesc(t,f,res.gPDC{mT,mD})
        set(gca,'YDir','normal')
        title(sprintf('%s_G_P_D_C %s -> %s',meth,siglab{mD},siglab{mT}));
        caxis([0 1])
        if(mD==1)
            ylabel('Hz')
        else
           set(gca, 'YTickLabelMode', 'Manual')
           set(gca, 'YTick', [])
        end
        if(mT==M)
            xlabel('n(sec)')
        else
           set(gca, 'XTickLabelMode', 'Manual')
           set(gca, 'XTick', [])
        end
        set(gca,'fontsize',12,'fontweight','bold')
    end
end    
saveas(gcf,sprintf('%sgPDC.png',meth))
