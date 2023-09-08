function [figHandle,figHandle2,figHandle3,figHandle4] = phiAndX(tPowerMin,tPowerMax,lPowerMin,lPowerMax,behaviorList,baseParams,nSamps)
    
    %Nicholas Szczecinski
    %West Virginia University
    %13 October 2020

    lPower = linspace(lPowerMin,lPowerMax,nSamps);
    tPower = linspace(tPowerMin,tPowerMax,nSamps);

    %Generate vectors for t and l
    l = 10.^lPower;
    t = 10.^tPower;

    %Generate 2D grids for t and l
    [T,L] = meshgrid(t,l);
    
    %Read in the parameter values when L=1.
    m0 = baseParams(1);
    c0 = baseParams(2);
    k0 = baseParams(3);
    s = baseParams(4);
    g = 10;
    
    %Because the torsional analogs of m, c, and k are complicated, define
    %them as their own function, which we can use consistently everywhere.
    elasOfL = @(L) k0*s^2*L.^3;
    gravOfL = @(L) m0*g/2*L.^4;
    mOfL = @(L) 1/3*m0*L.^5;
    cOfL = @(L) c0*s^2*L.^3;
    kOfL = @(L) elasOfL(L) + gravOfL(L);
    
    %To get a basic sampling of m, c, and k over the sample space, just
    %plug in L and save the output.
    m = mOfL(L);
    c = cOfL(L);
    k = kOfL(L);
    
    %Forcing frequency of the limb given the cycle period T.
    omega = 2*pi./T;
    
    %phi is the phase angle between the forcing function and the appendage
    %displacement. 
    %Textbook calculation: phi = atan2d(2*zeta.*r,1-r.^2);
    %Formal calculation:
    phi = atan2d(c.*omega,k-m.*omega.^2);
    
    %Define a function that computes the damping ratio zeta as a function
    %of L
    zetaOfL = @(L) cOfL(L)./(2*sqrt(kOfL(L).*mOfL(L)));
    
    %Numerically compute the length scale at which zeta=1, signifying the
    %boundary between underdamping and overdamping (i.e. critically
    %damped).
    LcritDamped = fzero(@(L)1-zetaOfL(L),[1e-3,1]);

    phiLcritDamped = atan2d(cOfL(LcritDamped).*(2*pi./t),kOfL(LcritDamped)-mOfL(LcritDamped).*(2*pi./t).^2);
    
    
    %Plot the figure of animal behaviors on top of the dynamic regimes.
    figHandle3 = figure;
    
    contour(T,L,phi,5:10:175)
    colormap autumn
    
    surfAx = gca;
    xlabel('Time-scale (s)')
    ylabel('Length-scale (m)')
    
    surfAx.XScale = 'log';
    surfAx.YScale = 'log';
    hold on

    xlim([10^tPowerMin,10^tPowerMax])
    ylim([10^lPowerMin,10^lPowerMax])
    view(0,90)
    
    cbar = colorbar;
    cbar.Limits = [0,180];
    cbar.Ticks = [0,90,180];
    set(cbar,'YTickLabel', []);
    hYLabel = ylabel(cbar, 'kinetic                viscous                quasi-static');     
    set(hYLabel,'Rotation',-90);
    hYLabel.Position(1) = 1.8;
    
    for i=1:size(behaviorList,1)
        plot3(behaviorList{i,2},[1,1]*behaviorList{i,3},[200,200],'k','linewidth',1.5)
        text(behaviorList{i,2}(2)*behaviorList{i,4},behaviorList{i,3}*behaviorList{i,5},200,behaviorList{i},'fontweight','bold','fontsize',8,'Color','black','fontname','Arial')
    end
    set(figHandle3,'Position',[480.2000e+000 153.8000e+000 517.7500  403.3000])

    text(10^1.3,10^-3.8,200,'quasi-static','fontname','Arial','fontweight','bold','fontsize',10,'color','black','HorizontalAlignment','right')
    text(10^-2.5,10^-3.8,200,'viscous','fontname','Arial','fontweight','bold','fontsize',10,'color','black','HorizontalAlignment','left')
    text(10^-2.5,10^0.8,200,'kinetic','fontname','Arial','fontweight','bold','fontsize',10,'color','black','HorizontalAlignment','left')
    tt = text(10^1.8,10^0.8,200,'underdamped','fontname','Arial','fontweight','bold','fontsize',10,'color','black','HorizontalAlignment','left');
    set(tt,'Rotation',-90)
    tt = text(10^1.8,10^-2.3,200,'overdamped','fontname','Arial','fontweight','bold','fontsize',10,'color','black','HorizontalAlignment','left');
    set(tt,'Rotation',-90)
    
    drawnow
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figHandle4 = figure;

    fadedAutumn = autumn;
    desiredAlpha = 0.6;
    fadedAutumn(:,2) = desiredAlpha*fadedAutumn(:,2) + (1 - desiredAlpha);
    fadedAutumn(:,3) = 1 - desiredAlpha;

%     sp1 = subplot(1,2,1);
%     sf1 = surf(T,L,phi,'edgecolor','none','facecolor','interp');
% 
%     surfAx = gca;
%     xlabel('T (s)')
%     ylabel('L (m)')
%     zlabel('\phi (deg)')
%     
%     surfAx.XScale = 'log';
%     surfAx.YScale = 'log';
%     surfAx.XMinorTick = 'off';
%     surfAx.YMinorTick = 'off';
%     surfAx.XMinorGrid = 'off';
%     surfAx.YMinorGrid = 'off';
%     hold on
%     
%     xlim([10^tPowerMin,10^tPowerMax])
%     ylim([10^lPowerMin,10^lPowerMax])
%     zlim([0,180])
%     
%     xticks(logspace(-3,2,6))
%     yticks(logspace(-5,1,7))
%     zticks(0:45:180)
% 
%     zticklabels({'0','45','90','135','180'})
%     
%     view(-20,60)
% 
%     colormap(surfAx,fadedAutumn)

    spContours = subplot(2,2,1);
    levs = 5:10:175;
    contour3(T,L,phi,levs,'linewidth',1);
    
    contAx = gca;
    
    contAx.XScale = 'log';
    contAx.YScale = 'log';
    grid off
    box off
    hold on
    
    xlim([10^tPowerMin,10^tPowerMax])
    ylim([10^lPowerMin,10^lPowerMax])
    zlim([0,180])
    
    xticks([])
    yticks([])
    zticks([])
    
    view(-20,60)
    
    colormap(contAx,autumn);
    
    drawnow

    cutthroughContours = subplot(2,2,2);

    plot3(t,LcritDamped+zeros(size(t)),phiLcritDamped,'k','linewidth',1)
    
    contAx = gca;
    
    contAx.XScale = 'log';
    contAx.YScale = 'log';
    grid off
    box off
    hold on
    
    xlim([10^tPowerMin,10^tPowerMax])
    ylim([10^lPowerMin,10^lPowerMax])
    zlim([0,180])
    
    xticks([])
    yticks([])
    zticks([])
    
    view(-20,60)

    cutthroughSurf = subplot(2,2,3);

    surf(T,L,phi,'EdgeColor','none','FaceColor','interp')
    colormap autumn

    contAx = gca;
    
    contAx.XScale = 'log';
    contAx.YScale = 'log';
    grid off
    box off
    hold on
    
    xlim([10^tPowerMin,10^tPowerMax])
    ylim([10^lPowerMin,10^lPowerMax])
    zlim([0,180])
    
    xticks([])
    yticks([])
    zticks([])
    
    view(-20,60)

    axisLabels = subplot(2,2,4);

    contAx = gca;
    
    contAx.XScale = 'log';
    contAx.YScale = 'log';
    contAx.XMinorTick = 'off';
    contAx.XMinorGrid = 'off';
    contAx.YMinorTick = 'off';
    contAx.YMinorGrid = 'off';
    contAx.ZMinorTick = 'off';
    contAx.ZMinorGrid = 'off';
    grid on
    hold on
    
    xlim([10^tPowerMin,10^tPowerMax])
    ylim([10^lPowerMin,10^lPowerMax])
    zlim([0,180])
    
    tPowerMinRnd = ceil(tPowerMin);
    tPowerMaxRnd = ceil(tPowerMax);
    lPowerMinRnd = ceil(lPowerMin);
    lPowerMaxRnd = floor(lPowerMax);

    xticks(logspace(tPowerMinRnd,tPowerMaxRnd,tPowerMaxRnd - tPowerMinRnd + 1))
    yticks(logspace(lPowerMinRnd,lPowerMaxRnd,lPowerMaxRnd - lPowerMinRnd + 1))
    zticks(0:45:180)

    xlabel('T (s)')
    ylabel('L (m)')
    zlabel('\phi (deg)')
    
    view(-20,60)
    
    figHandle4.Position = [281 57 766 705];
    
    drawnow
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Define time and length scales of interest, which we will label in our
    %T and L plot with the letters A, B, C, D, and E.
    Lletter(1) = 3e-3;
    Tletter(1) = 3;
    
    Lletter(2) = 1.7;
    Tletter(2) = 3;
    
    Lletter(3) = 1.7;
    Tletter(3) = 1.5; %0.5;

    Lletter(4) = 2e-3;
    Tletter(4) = 1e-2;

    Lletter(5) = .03;
    Tletter(5) = 10^-2.5;
    
    letterCell = {'D','C','B','E','F'};
    letterColor = {'black','black','black','black','black'};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figHandle = figure;
    spSurf = subplot(3,3,[1,2,4,5]);
    contour(T,L,phi,5:10:175)
    colormap autumn
    surfAx = gca;
    xlabel('Time-scale (s)')
    ylabel('Length-scale (m)')
    title('A')
    
    surfAx.XScale = 'log';
    surfAx.YScale = 'log';
    % surfAx.ZScale = 'log';
    hold on
    lw = 1;
    bndColor = 0*[1,1,1];
%     plot3(TnatFreq,LnatFreq,180+zeros(size(TnatFreq)),'color',bndColor,'linewidth',lw,'linestyle','-')
%     plot3(Tev,Lev,180+zeros(size(Tev)),'color',bndColor,'linewidth',lw,'linestyle','--')
%     plot3(Tvi,Lvi,180+zeros(size(Tvi)),'color',bndColor,'linewidth',lw,'linestyle','-.')
    plot3([10^tPowerMin,10^tPowerMax],LcritDamped*[1,1],[180,180],'color',bndColor,'linewidth',lw,'linestyle',':')
    xlim([10^tPowerMin,10^tPowerMax])
    ylim([10^lPowerMin,10^lPowerMax])
    view(0,90)
    
    cbar = colorbar;
    cbar.Limits = [0,180];
    cbar.Ticks = [0,90,180];
    set(cbar,'YTickLabel', []);
    hYLabel = ylabel(cbar, 'kinetic          viscous           quasi-static');     
    set(hYLabel,'Rotation',-90);
    hYLabel.Position(1) = 1.8;
    cbar.Position(1) = .53;

    set(figHandle,'Position',[480.2000e+000 153.8000e+000 517.7500  403.3000])
    
    surfAx.XTick = logspace(tPower(1),tPower(end),tPower(end)-tPower(1)+1);
    surfAx.YTick = logspace(lPower(1),lPower(end),lPower(end)-lPower(1)+1);
    surfAx.Position(4) = surfAx.Position(4) - 0.05;
    surfAx.Position(2) = surfAx.Position(2) + 0.05;
    cbar.Position(2) = surfAx.Position(2);
    cbar.Position(4) = surfAx.Position(4);
    
    text(10^1.3,10^-3.8,200,'quasi-static','fontname','Arial','fontweight','bold','fontsize',8,'color','black','HorizontalAlignment','right')
    text(10^-2.5,10^-3.8,200,'viscous','fontname','Arial','fontweight','bold','fontsize',8,'color','black','HorizontalAlignment','left')
    text(10^-2.5,10^0.8,200,'kinetic','fontname','Arial','fontweight','bold','fontsize',8,'color','black','HorizontalAlignment','left')
    tt = text(10^1.8,10^0.8,200,'underdamped','fontname','Arial','fontweight','bold','fontsize',8,'color','black','HorizontalAlignment','left');
    set(tt,'Rotation',-90)
    tt = text(10^1.8,10^-1.8,200,'overdamped','fontname','Arial','fontweight','bold','fontsize',8,'color','black','HorizontalAlignment','left');
    set(tt,'Rotation',-90)
    
    
    sp1 = [];
    sp2 = [];

    tSim = cell(5,1);
    xSim = cell(5,1);
    Fapp = cell(5,1);

    for i=1:5
        [tSim{i},xSim{i},~,~,~,~,~,Fapp{i},hFig] = simulateJointResponse(k0,c0,m0,s,Lletter(i),3,Tletter(i),15*Tletter(i),sp1,sp2);
        close(hFig)
    end
    
    spInd = [9,6,3,8,7];
    
    figure(figHandle)
    
    for i=1:5
        %Label on the T vs. L plot
        subplot(spSurf)
        hold on
        box on
        text(Tletter(i),Lletter(i),200,letterCell{i},'color',letterColor{i},'HorizontalAlignment','center','VerticalAlignment','middle','fontweight','bold')
        
        %Pick the subplot
        spWL = subplot(3,3,spInd(i));
        
        %Plot work loop
        title(letterCell{i})
        pts = (tSim{i} <= Tletter(i) & Fapp{i}(tSim{i}) > 0);
        hold on
        
        cvec = autumn(3);
        %x sim is negative because of how a muscle would shorten. Consider:
        %extensor causes positive torque. But when the joint is moving in
        %the positive direction, the muscle is shortening (negative
        %direction). Thus the joint angle should have the opposite sign to
        %measure the work done by the actuator.   
        
        time = tSim{i}(pts);
        F = max(Fapp{i}(time),0);
        x = -xSim{i}(pts);
        
        [maxF,ind] = max(F);
        if x(ind) < -0.1
            %quasi-static
        
            %Find the index at which the minimum position is reached. Split
            %up F and x before and after this point.
            [~,ind] = min(x);
            ftop = F(1:ind);
            xtop = x(1:ind);
            fbot = F(ind:end);
            xbot = x(ind:end);

            %Sample the "top" and "bottom" of the loop at the same grid,
            %allowing us to use the area() plot.
            Xbot = -.5:.01:.5;
            Fbot = interp1(xbot,fbot,Xbot);
            Fbot(isnan(Fbot)) = 0;
            Xtop = Xbot;
            Ftop = interp1(xtop,ftop,Xtop);
            Ftop(isnan(Ftop)) = 0;

            %The extremum can't be mapped with interp1. Manually set these
            %values.
            Fbot(Xbot == -1) = max(F);
            Ftop(Xtop == -1) = max(F);
            
            %Plot the area under these curves and color code them based on
            %where the energy goes.
            area(Xbot',Ftop','edgealpha',0,'facecolor',cvec(2,:));
            area(Xbot',Fbot','edgealpha',0,'facecolor',cvec(1,:));
            plot(Xbot,Fbot,'k','linewidth',1)
            plot(Xtop,Ftop,'k','linewidth',1)
            
            %Draw directional arrows
            botArrLoc = ceil(length(xbot)/2);
            botStrt = [xbot(botArrLoc),fbot(botArrLoc)];
            botStop = [xbot(botArrLoc+1),fbot(botArrLoc+1)];
            arrObj = arrow(botStrt,botStop);
            
            topArrLoc = ceil(length(xtop)/3);
            botStrt = [xtop(topArrLoc),ftop(topArrLoc)];
            botStop = [xtop(topArrLoc+1),ftop(topArrLoc+1)];
            arrObj = arrow(botStrt,botStop);
        elseif x(ind) > 0.1
            %inertial (kinetic)
            
            %Find the index at which the minimum position is reached. Split
            %up F and x before and after this point.
            [~,ind] = max(x);
            ftop = F(ind:end);
            xtop = x(ind:end);
            fbot = F(1:ind);
            xbot = x(1:ind);

            %Sample the "top" and "bottom" of the loop at the same grid,
            %allowing us to use the area() plot.
            Xbot = -.5:.01:.5;
            Fbot = interp1(xbot,fbot,Xbot);
            Fbot(isnan(Fbot)) = 0;
            Xtop = Xbot;
            Ftop = interp1(xtop,ftop,Xtop);
            Ftop(isnan(Ftop)) = 0;
            
            %The extremum can't be mapped with interp1. Manually set these
            %values.
            Fbot(Xbot == 1) = max(F);
            Ftop(Xtop == 1) = max(F);

            %Plot the area under these curves and color code them based on
            %where the energy goes.
            area(Xbot',Ftop','edgealpha',0,'facecolor',cvec(2,:));
            area(Xbot',Fbot','edgealpha',0,'facecolor',cvec(3,:));
            plot(Xbot',Ftop','k','linewidth',1)
            plot(Xbot',Fbot','k','linewidth',1)
            
            %Draw directional arrows
            botArrLoc = ceil(length(xbot)/3);
            botStrt = [xbot(botArrLoc),fbot(botArrLoc)];
            botStop = [xbot(botArrLoc+1),fbot(botArrLoc+1)];
            arrObj = arrow(botStrt,botStop);
            
            topArrLoc = ceil(length(xtop)/2);
            botStrt = [xtop(topArrLoc),ftop(topArrLoc)];
            botStop = [xtop(topArrLoc+1),ftop(topArrLoc+1)];
            arrObj = arrow(botStrt,botStop);
        else
            %viscous
            
            [~,ind] = min(x);
            ftop = F(1:ind);
            xtop = x(1:ind);
            fbot = F(ind:end);
            xbot = x(ind:end);

            Xbot = -.5:.01:.5;
            Fbot = zeros(size(Xbot));
            Xtop = Xbot;
            Ftop = interp1(xtop,ftop,Xtop);
            Ftop(isnan(Xtop)) = 0;
            
            %The extremum can't be mapped with interp1. Manually set these
            %values.
            Fbot(Xbot == 1) = 0;
            Ftop(Xtop == 1) = 0;
            Fbot(Xbot == -1) = 0;
            Ftop(Xbot == -1) = 0;

            area(Xtop',Ftop','edgealpha',0,'facecolor',cvec(2,:));
            plot(Xtop',Ftop','k','linewidth',1)
            plot(Xbot',Fbot','k','linewidth',1)
            
            %Draw directional arrows
            botArrLoc = ceil(length(xbot)/2);
            botStrt = [-1,0];
            botStop = [.1,0];
            arrObj = arrow(botStrt,botStop);
            arrow fixlimits
            
            topArrLoc = ceil(length(xtop)/2);
            botStrt = [xtop(topArrLoc),ftop(topArrLoc)];
            botStop = [xtop(topArrLoc+1),ftop(topArrLoc+1)];
            arrObj = arrow(botStrt,botStop);
        end        

        if ~strcmp(letterCell{i},'B') && ~strcmp(letterCell{i},'C')
            xlabel('Angle, \theta (rad)')
        end
        if ~strcmp(letterCell{i},'D') && ~strcmp(letterCell{i},'E')
            ylabel('Moment, M (Nm)')
        end
        
        title(letterCell{i})
        
        box on
        set(gca, 'Layer', 'Top')
        
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figHandle2 = figure;
    sp = subplot(3,3,[1,2,4,5]);
    contour(T,L,phi,5:10:175)
    colormap autumn
    surfAx = gca;
    xlabel('Time-scale (s)')
    ylabel('Length-scale (m)')
    title('A')
    
    surfAx.XScale = 'log';
    surfAx.YScale = 'log';
    hold on
    lw = 1;
    bndColor = 0*[1,1,1];
   

    plot3([10^tPowerMin,10^tPowerMax],LcritDamped*[1,1],[180,180],'color',bndColor,'linewidth',lw,'linestyle',':')
    xlim([10^tPowerMin,10^tPowerMax])
    ylim([10^lPowerMin,10^lPowerMax])
    view(0,90)
    
    cbar = colorbar;
    cbar.Limits = [0,180];
    cbar.Ticks = [0,90,180];
    set(cbar,'YTickLabel', []);
    hYLabel = ylabel(cbar, 'kinetic              viscous       quasi-static');     
    set(hYLabel,'Rotation',-90);
    hYLabel.Position(1) = 1.8;
    cbar.Position(1) = .53;
    

    set(figHandle2,'Position',[480.2000e+000 153.8000e+000 517.7500  403.3000])
    surfAx.Position(4) = surfAx.Position(4) - 0.05;
    surfAx.Position(2) = surfAx.Position(2) + 0.05;
    cbar.Position(2) = surfAx.Position(2);
    cbar.Position(4) = surfAx.Position(4);
    
    surfAx.XTick = logspace(tPower(1),tPower(end),tPower(end)-tPower(1)+1);
    surfAx.YTick = logspace(lPower(1),lPower(end),lPower(end)-lPower(1)+1);

    text(10^1.3,10^-3.8,200,'quasi-static','fontname','Arial','fontweight','bold','fontsize',8,'color','black','HorizontalAlignment','right')
    text(10^-2.5,10^-3.8,200,'viscous','fontname','Arial','fontweight','bold','fontsize',8,'color','black','HorizontalAlignment','left')
    text(10^-2.5,10^0.8,200,'kinetic','fontname','Arial','fontweight','bold','fontsize',8,'color','black','HorizontalAlignment','left')
    tt = text(10^1.8,10^0.8,200,'underdamped','fontname','Arial','fontweight','bold','fontsize',8,'color','black','HorizontalAlignment','left');
    set(tt,'Rotation',-90)
    tt = text(10^1.8,10^-1.8,200,'overdamped','fontname','Arial','fontweight','bold','fontsize',8,'color','black','HorizontalAlignment','left');
    set(tt,'Rotation',-90)
    
    stimMultiplier = [1;1;3;1;1]+.5;
    numCycles = [5;10;10;10;20];
    
    sp2 = [];
    
    for i=1:5
        sp1 = subplot(3,3,spInd(i));
        [tSim{i},xSim{i},~,~,~,~,~,Fapp{i},~] = simulateJointResponse(k0,c0,m0,s,Lletter(i),numCycles(i),Tletter(i),stimMultiplier(i)*Tletter(i),sp1,sp2);
        title(letterCell{i})
        
        if ~strcmp(letterCell{i},'B') && ~strcmp(letterCell{i},'C')
            xlabel('Time, t (n.d.)')
        else
            xlabel('')
        end
        if ~strcmp(letterCell{i},'D') && ~strcmp(letterCell{i},'E')
            ylabel('Angle, \theta (rad)')
        else
            ylabel('')
        end
    end
    
    %Add titles to subplot
    subplot(sp)
    for i=1:5
        text(Tletter(i),Lletter(i),200,letterCell{i},'color',letterColor{i},'HorizontalAlignment','center','VerticalAlignment','middle','fontweight','bold')
    end   

    hCutthroughSummary = figure;
    hCutthroughSummary.Position = [488 41.8000 315.4000 740.8000];

    subplot(6,1,3:4)

    surf(T,L,phi,'edgecolor','none','facecolor','interp')
    hold on
%     contour3(T,L,phi,5:10:175,'linewidth',1)
%     ch.ZLocation = -90;
    colormap(fadedAutumn)
    
    surfAx = gca;
    xlabel('T (s)')
    ylabel('L (m)')
    zlabel('\phi (deg)')
    
    surfAx.XScale = 'log';
    surfAx.YScale = 'log';
    surfAx.XMinorGrid = 'off';
    surfAx.YMinorGrid = 'off';
    surfAx.XMinorTick = 'off';
    surfAx.YMinorTick = 'off';
    hold on
    
    xlim([10^tPowerMin,10^tPowerMax])
    ylim([10^lPowerMin,10^lPowerMax])
    zlim([0,180])
    
    xticks(logspace(-3,2,6))
    yticks(logspace(-5,1,7))
    zticks(0:45:180)
    zticklabels({'0','45','90','135','180'})
    
    view(-20,60)

    LindToPlot = 1;
    TindToPlot = 2;
    
    %Calculate relative magnitude of forces at constant length as frequency
    %changes
    Lscale = [5e-3,3e-2,1]; %m
    om = 2*pi./t;
    for i=1:length(Lscale)
        torqueInertialLength = mOfL(Lscale(i))*om.^2;
        torqueViscousLength = cOfL(Lscale(i))*om;
        torqueElasticLength = elasOfL(Lscale(i)) + zeros(size(om));
        torqueGravLength = gravOfL(Lscale(i)) + zeros(size(om));
    
        torqueNormLength = sqrt(torqueInertialLength.^2 + torqueViscousLength.^2 + torqueElasticLength.^2 + torqueGravLength.^2);
    
        m = mOfL(Lscale(i));
        c = cOfL(Lscale(i));
        k = kOfL(Lscale(i));
        
        %phi is the phase angle between the forcing function and the appendage
        %displacement. 
        %Textbook calculation: phi = atan2d(2*zeta.*r,1-r.^2);
        %Formal calculation:
        phiLength = atan2d(c.*om,k-m.*om.^2);
    
        figure
        subplot(2,1,1)
        plot(t,torqueElasticLength./torqueNormLength,'linewidth',1,'color',cvec(1,:));
        hold on
        plot(t,torqueGravLength./torqueNormLength,'--','linewidth',1,'color',cvec(1,:));
        plot(t,torqueViscousLength./torqueNormLength,'linewidth',1,'color',cvec(2,:));
        plot(t,torqueInertialLength./torqueNormLength,'linewidth',1,'color',cvec(3,:));
        ax = gca;
        ax.XScale = 'log';
        ylabel('relative magnitude')
        grid on
        legend('elas.','grav.','visc.','iner.')
        
        subplot(2,1,2)
        plot(t,phiLength,'k','linewidth',1)
        ax = gca;
        ax.XScale = 'log';
        xlabel('T (s)')
        ylabel('\phi (deg)')
        ylim([0,180]);
        yticks(0:90:180)
        grid on
    
        sgtitle(sprintf('Relative force magnitude and phase shift in swing, L = %0.3f m',Lscale(i)))

        if i == LindToPlot
            figure(hCutthroughSummary)

            sp5 = subplot(6,1,5);
            plot(t,torqueElasticLength./torqueNormLength,'linewidth',1,'color',cvec(1,:));
            hold on
            plot(t,torqueGravLength./torqueNormLength,'--','linewidth',1,'color',cvec(1,:));
            plot(t,torqueViscousLength./torqueNormLength,'linewidth',1,'color',cvec(2,:));
            plot(t,torqueInertialLength./torqueNormLength,'linewidth',1,'color',cvec(3,:));
            ax = gca;
            ax.XScale = 'log';
            ylabel('relative magnitude')
            grid on
            legend('elas.','grav.','visc.','iner.')
            xlim([min(t),max(t)])
            ax.XMinorGrid = 'off';
            ax.YMinorGrid = 'off';
            ax.XMinorTick = 'off';
            ax.YMinorTick = 'off';
            xticks(logspace(log10(min(t)),log10(max(t)),log10(max(t))-log10(min(t))+1))
            set(sp5,'Position',[0.1300 0.2525 0.7750 0.1026]);
            title(sprintf('C. Relative force magnitude and phase shift in swing, L = %0.3f m',Lscale(i)),'FontSize',8,'FontWeight','normal')
            
            sp6 = subplot(6,1,6);
            plot(t,phiLength,'k','linewidth',1)
            ax = gca;
            ax.XScale = 'log';
            xlabel('T (s)')
            ylabel('\phi (deg)')
            ylim([0,180]);
            yticks(0:90:180)
            ax.XMinorGrid = 'off';
            ax.YMinorGrid = 'off';
            ax.XMinorTick = 'off';
            ax.YMinorTick = 'off';
            grid on
            xlim([min(t),max(t)])
            xticks(logspace(log10(min(t)),log10(max(t)),log10(max(t))-log10(min(t))+1))
            set(sp6,'Position',[0.1300 0.1100 0.7750 0.1026]);

            subplot(6,1,3:4)
            hold on
            plot3(t,Lscale(LindToPlot)+zeros(size(t)),phiLength+.1,'k-','linewidth',2)
        
%             figure(figHandle4)
%             subplot(cutthroughContours)
%             hold on
%             plot3(t,Lscale(LindToPlot)+zeros(size(t)),phiLength+.1,'k-','linewidth',2)
        
        end
    end
    %Calculate relative magnitude of forces at constant frequency as length 
    % scale changes
    Tscale = [0.1,1,10]; %s
    om = 2*pi./Tscale;
    for i=1:length(Tscale)
        torqueInertialFreq = mOfL(l)*om(i).^2;
        torqueViscousFreq = cOfL(l)*om(i);
        torqueElasticFreq = elasOfL(l);
        torqueGravFreq = gravOfL(l);
    
        torqueNormFreq = sqrt(torqueInertialFreq.^2 + torqueViscousFreq.^2 + torqueElasticFreq.^2 + torqueGravFreq.^2);
    
        m = mOfL(l);
        c = cOfL(l);
        k = kOfL(l);
        
        %phi is the phase angle between the forcing function and the appendage
        %displacement. 
        %Textbook calculation: phi = atan2d(2*zeta.*r,1-r.^2);
        %Formal calculation:
        phiFreq = atan2d(c.*om(i),k-m.*om(i).^2);
    
        figure
        subplot(2,1,1)
        plot(l,torqueElasticFreq./torqueNormFreq,'linewidth',1,'color',cvec(1,:));
        hold on
        plot(l,torqueGravFreq./torqueNormFreq,'--','linewidth',1,'color',cvec(1,:));
        plot(l,torqueViscousFreq./torqueNormFreq,'linewidth',1,'color',cvec(2,:));
        plot(l,torqueInertialFreq./torqueNormFreq,'linewidth',1,'color',cvec(3,:));
        ax = gca;
        ax.XScale = 'log';
        xlim([min(l),max(l)])
        ylabel('relative magnitude')
        grid on
        legend('elas.','grav.','visc.','iner.')
        
        subplot(2,1,2)
        plot(l,phiFreq,'k','linewidth',1)
        ax = gca;
        ax.XScale = 'log';
        xlabel('L (m)')
        ylabel('\phi (deg)')
        ylim([0,180]);
        yticks(0:90:180)
        xlim([min(l),max(l)])
        grid on
    
        sgtitle(sprintf('Relative force magnitude and phase shift in swing, T = %0.3f s',Tscale(i)))

        if i == TindToPlot
            figure(hCutthroughSummary)

            sp1 = subplot(6,1,1);
            plot(l,torqueElasticFreq./torqueNormFreq,'linewidth',1,'color',cvec(1,:));
            hold on
            plot(l,torqueGravFreq./torqueNormFreq,'--','linewidth',1,'color',cvec(1,:));
            plot(l,torqueViscousFreq./torqueNormFreq,'linewidth',1,'color',cvec(2,:));
            plot(l,torqueInertialFreq./torqueNormFreq,'linewidth',1,'color',cvec(3,:));
            ax = gca;
            ax.XScale = 'log';
            xlim([min(l),max(l)])
            xticks(logspace(ceil(log10(min(l))),floor(log10(max(l))),floor(log10(max(l)))-ceil(log10(min(l)))+1))
            ylabel('relative magnitude')
            ax.XMinorGrid = 'off';
            ax.YMinorGrid = 'off';
            ax.XMinorTick = 'off';
            ax.YMinorTick = 'off';
            grid on
            legend('elas.','grav.','visc.','iner.')
            title(sprintf('A. Relative force magnitude and phase shift in swing, T = %0.1f s',Tscale(i)),'FontSize',8,'FontWeight','normal')
            
            
            sp2 = subplot(6,1,2);
            plot(l,phiFreq,'k:','linewidth',1.5)
            ax = gca;
            ax.XScale = 'log';
            xlabel('L (m)')
            ylabel('\phi (deg)')
            ylim([0,180]);
            yticks(0:90:180)
            xlim([min(l),max(l)])
            xticks(logspace(ceil(log10(min(l))),floor(log10(max(l))),floor(log10(max(l)))-ceil(log10(min(l)))+1))
            ax.XMinorGrid = 'off';
            ax.YMinorGrid = 'off';
            ax.XMinorTick = 'off';
            ax.YMinorTick = 'off';
            grid on
            

            subplot(6,1,3:4)
            hold on
            plot3(Tscale(TindToPlot)+zeros(size(l)),l,phiFreq+.1,'k:','linewidth',2)
            zlim([0,180])
            zticks(0:45:180)

%             figure(figHandle4)
%             subplot(cutthroughContours)
%             hold on
%             plot3(Tscale(TindToPlot)+zeros(size(l)),l,phiFreq+.1,'k:','linewidth',2)
            
        end
    end

    set(figHandle4,'renderer','painters')
    set(hCutthroughSummary,'renderer','painters') 
    

end