function [t,x,v,Fsp,Fg,Fvis,Finer,Fapp,h,TapplyPert] = simulateJointResponseStance(k0,c0,m0,s,a,L,numCycles,T,Tperturb,spObj1,spObj2)
    
    %Nicholas Szczecinski
    %West Virginia University
    %13 October 2020

    g = 10;
    kSp = k0*s^2*L.^3;
    kG = -(a + 1/2)*m0*g*L.^4;
    k = kSp + kG;
    J = (a + 1/3)*m0*L.^5;
    c = c0*s^2*L.^3;
    tmax = numCycles*T;

    omega = 2*pi/T;
    
    X = 0.5; %1; %We are enforcing this be true, based on our definition of F0 below.
    phi = atan(c*omega/(k - J*omega^2));
    while phi < 0
        phi = phi + pi;
    end
    while phi > pi
        phi = phi - pi;
    end
    
    xp0 = X*sin(-phi);
    xp_dot0 = omega*X*cos(-phi);
    
    F0 = X*sqrt((k - J*omega^2)^2 + c^2*omega^2);
    
    %To make sure perturbation is always applied at the same phase, add the
    %phase lag to the commanded perturbation time
    TapplyPert = Tperturb - T*phi/(2*pi);
    pertDuration = T/2;
    
    M = [1,0;0,J];
        
    options = odeset('Mass',M,'RelTol',100*eps,'AbsTol',eps,'Events',@(t,y)fQuit(t,y));
    y0 = [xp0;xp_dot0];
%     tspan = [0,tmax];
    tspan = 0:T/1000:tmax;
    if tspan(end) ~= tmax
        tspan(end+1) = tmax;
    end
    
    Fapp = @(t)F0*sin(omega*t);
    f = @(t,x) [x(2);-k*x(1) - c*x(2) + Fapp(t)];
    [t1,y1,te1,ye1,ie1] = ode15s(f,[0,TapplyPert],y0,options);
    
    if isempty(ie1) && t1(end) < tmax
        %The simulation completed without reaching an event, so we can
        %continue it with the perturbation force.
        fprintf('start pert\n')
        y0 = y1(end,:)';
        Fapp = @(t)F0*sin(omega*t) + F0/5;
        f = @(t,x) [x(2);-k*x(1) - c*x(2) + Fapp(t)];
        [t2,y2,te2,ye2,ie2] = ode15s(f,[TapplyPert,TapplyPert + pertDuration],y0,options);
        
        if isempty(ie2) && t2(end) < tmax
            %The simulation completed without reaching an event, so we can
            %continue it without the perturbation force.
            y0 = y2(end,:)';
            Fapp = @(t)F0*sin(omega*t);
            f = @(t,x) [x(2);-k*x(1) - c*x(2) + Fapp(t)];
            if tmax > TapplyPert + pertDuration
                tspan = [TapplyPert + pertDuration,tmax];
            else
                tspan = [TapplyPert + pertDuration, TapplyPert + 3*T];
            end
            [t3,y3,te3,ye3,ie3] = ode15s(f,tspan,y0,options);

            x = [y1(:,1);y2(:,1);y3(:,1)];
            v = [y1(:,2);y2(:,2);y3(:,2)];
            t = [t1;t2;t3];
        else
            x = [y1(:,1);y2(:,1)];
            v = [y1(:,2);y2(:,2)];
            t = [t1;t2];
        end
    else
        %The simulation completed upon reaching an event, so the simulation
        %should be stopped.
        x = y1(:,1);
        v = y1(:,2);
        t = t1;
    end
    
    a = centeredDiff(t,v);
    
    Fsp = -kSp*x;
    Fg = -kG*x;
    Fvis = -c*v;
    Finer = -J*a;
    
    if ~isempty(spObj1) && ~isempty(spObj2)
        h = get(spObj1,'parent');
        hold on
        newFig = false;
    elseif isempty(spObj1) && isempty(spObj2)
        h = figure;
        spObj1 = subplot(2,1,1);
        spObj2 = subplot(2,1,2);
        newFig = true;
    elseif ~isempty(spObj1) && isempty(spObj2)
        h = get(spObj1,'parent');
        spObj2 = [];
        newFig = false;
    elseif isempty(spObj1) && ~isempty(spObj2)
        h = get(spObj2,'parent');
        spObj1 = [];
        newFig = false;
    end
    
    if ~isempty(spObj1)
        subplot(spObj1)
        plot(t,x,'k','linewidth',1)
        if isempty(spObj2)
            xlabel('t')
            if numCycles <= 5
                xticks(0:T:tmax)
                xticklabels({'0','T','2T','3T','4T','5T'})
            elseif numCycles <= 10
                xticks(0:2*T:tmax)
                xticklabels({'0','2T','4T','6T','8T','10T'})
            else
                xticks(0:4*T:tmax)
                xticklabels({'0','4T','8T','12T','16T','20T'})
            end
            xlim([0,tspan(end)])
        else
            xticks([])
        end
        ylabel('x')
        
        if TapplyPert < tmax
            hold on
            plot([TapplyPert,min(tmax,TapplyPert + pertDuration)],[-X-.5,-X-.5],'k','linewidth',2)
        end
        
    end
    
    if ~isempty(spObj2)
        subplot(spObj2)
%         colorVec = jet(5);
        colorVec = parula(5);
        plot(t,Fapp(t),'k','linewidth',1)
        hold on
        if newFig
            plot(t,Fsp,'color',colorVec(1,:),'linewidth',1);
            plot(t,Fg,'--','color',colorVec(1,:),'linewidth',1);
            plot(t,Fvis,'color',colorVec(3,:),'linewidth',1);
            plot(t,Finer,'color',colorVec(5,:),'linewidth',1);
            legend('applied','potential, elas','potential, grav','viscous','inertial')
        else
            plot(t,Fsp+Fg,'color',colorVec(1,:),'linewidth',1);
            plot(t,Fvis,'color',colorVec(3,:),'linewidth',1);
            plot(t,Finer,'color',colorVec(5,:),'linewidth',1);
        end


        xlabel('t')
        if numCycles <= 5
            xticks(0:T:tmax)
            xticklabels({'0','T','2T','3T','4T','5T'})
        elseif numCycles <= 10
            xticks(0:2*T:tmax)
            xticklabels({'0','2T','4T','6T','8T','10T'})
        else
            xticks(0:4*T:tmax)
            xticklabels({'0','4T','8T','12T','16T','20T'})
        end
        xlim([tspan(1),tspan(end)])
        ylabel('F')
    end
    
    drawnow
end