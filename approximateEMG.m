function [agEMG,antagEMG,tVec,stMask] = approximateEMG(tVec,theta,tStanceToSwing,omegaEMG,T,L,baseParams)
    %Read in the parameter values when L=1.
    m0 = baseParams(1);
    c0 = baseParams(2);
    k0 = baseParams(3);
    s = baseParams(4);
    a = baseParams(5);
    thEq = baseParams(6);
    g = 10;
    
    if iscolumn(theta)
        %good
    elseif isrow(theta)
        %transpose
        theta = theta';
    else
        error('theta must be a column vector.')
    end
    
    mSt = (a + 1/3)*m0*L.^5;
    cSt = c0*s^2*L.^3;
    kStElas = k0*s^2*L.^3;
    kStGrav = -(a + 1/2)*m0*g*L.^4;
    
    mSw = 1/3*m0*L.^5;
    cSw = c0*s^2*L.^3;
    kSwElas = k0*s^2*L.^3;
    kSwGrav = m0*g/2*L.^4;
    
    thDot = centeredDiff(tVec,theta);
    thDDot = centeredDiff(tVec,thDot);
    
    %We need swing and stance masks that repeat with each cycle.
    stMask = false(size(tVec));
    nCycles = ceil(max(tVec)/T);
    
    for i=1:nCycles
        stMask(tVec <= (i-1)*T + tStanceToSwing & tVec > (i-1)*T) = true;
    end
    swMask = ~stMask;
    
    FstanceInertia = mSt*thDDot(stMask);
    FswingInertia = mSw*thDDot(swMask);
    
    FstanceViscous = cSt*thDot(stMask);
    FswingViscous = cSw*thDot(swMask);
    
    FstanceElas = kStElas*(theta(stMask)-thEq);
    FswingElas = kSwElas*(theta(swMask)-thEq);
    
    FstanceGrav = kStGrav*theta(stMask);
    FswingGrav = kSwGrav*theta(swMask);
    
    Fst = FstanceInertia + FstanceViscous + FstanceElas + FstanceGrav;
    Fsw = FswingInertia + FswingViscous + FswingElas + FswingGrav;
    
    F = NaN(size(tVec));
    F(stMask) = Fst;
    F(swMask) = Fsw;
    
    
    Finertia(stMask) = FstanceInertia;
    Finertia(swMask) = FswingInertia;
    
    Fviscous(stMask) = FstanceViscous;
    Fviscous(swMask) = FswingViscous;
    
    Felastic(stMask) = FstanceElas;
    Felastic(swMask) = FswingElas;
    
    Fgrav(stMask) = FstanceGrav;
    Fgrav(swMask) = FswingGrav;
    
    figure
    plot(tVec,-Finertia,'linewidth',1)
    hold on
    plot(tVec,-Fviscous,'linewidth',1)
    plot(tVec,-Felastic,'linewidth',1)
    plot(tVec,-Fgrav,'linewidth',1)
    plot(tVec,F,'k','linewidth',1)
    legend('inertia','viscous','elastic','grav','net')
    xlabel('time (s)')
    ylabel('torque (Nm)')
    title('forces experienced during locomotion')
    xlim([0,max(tVec)])
    
    Fagonist = max(F,0);
    Fantagonist = -min(F,0);
    
    if ~isempty(omegaEMG)
        %If a value is provided for the cutoff frequency of the EMG low
        %pass filter, then filter the force.
        options = odeset('Mass',omegaEMG,'RelTol',100*eps,'AbsTol',eps);        
        
        f = @(t,x) interp1(tVec,Fagonist,t) - x(1);
        [tVec,agEMG] = ode45(f,tVec,Fagonist(1),options);
        
        f = @(t,x) interp1(tVec,Fantagonist,t) - x(1);
        [tVec,antagEMG] = ode45(f,tVec,Fantagonist(1),options);
    else
        %If no value is provided for the cutoff frequency of the EMG low
        %pass filter, then simply return the calculated force as the EMG
        %pattern.
        agEMG = Fagonist;
        antagEMG = Fantagonist;
    end
    
    figure
    plot(tVec,Fagonist)
    hold on
    plot(tVec,agEMG)

    figure
    plot(tVec,Fantagonist)
    hold on
    plot(tVec,antagEMG)
end