clear
close all


m0 = 12;
c0 = 1.31e3; %see zetaAcrossScales.m
k0 = 12e3; 
s = sqrt(1e-3);
a = 10;
lw = 1;

%Range of motion of the joint
h = 1; %rad

%Test sawtooth motion at different speeds for one size.
n = 1000;

omegaEMG = .02; %[]; %.02; %2*pi/30; %30 Hz filter means omega = 2*pi/f.
% in Hooper et al. 2009, tau = 0.02 s for the stick insect EMG.

tEMGdelay = 0.05; %0.05; %50 ms delay between EMG production and the generation of
%force. Because our model outputs force, we must shift each curve leftward
%by 50 ms relative to stance and swing phase.

%person
%Nilsson, Thorstensson, and Halbertsma 1985
% speeds = [0.6, 1.6, 3, 4, 8]; %m/s %walk walk walk run run
% % T = [1.8, 1, 0.75, .7, .55];
% % durSwing = [.5,.35,.3,.45,.4];
L = 1;
thEq = 0;

%speed: 3 m/s
T = 0.75;
durSwing = 0.3;

%speed: 8 m/s
% T = 0.55;
% durSwing = 0.4;

tStanceToSwing = T - durSwing;

%cat
%Smith, Chung, and Zernicke 1993
% speeds = [0.5,0.8,1.75]; %m/s
% T = [.7,.55,.35];
% T = .55;
% L = 15e-2;
% durSwing = 250e-3;
% tStanceToSwing = T - durSwing;

% %fruit fly:
%Szczecinski, Bockemuehl et al. 2018
% % T = linspace(0.06,.25,5);
T = .1;
L = 0.5e-3; 
durSwing = 0.03;
tStanceToSwing = T - durSwing;
tEMGdelay = 8e-3; %Azevedo, Dickinson, et al. 2020, fruit fly EMG
omegaEMG = .002;

% %cockroach:
% T = linspace(0.125,.375,3);
% T = .25;
% L = 8e-3; 
% tStanceToSwing = T - 1/2*min(T);

% %stick insect: Compare to Dallmann, Duerr, and Schmitz 2019
% %Figure 4B shows EMGs from protractor and retractor of freely walking stick insects
% T = 1;
% L = 15e-3; 
% tStanceToSwing = .75;
% thEq = -.25;

% %horse overground walk, Table 2 and Fig. 8 of Harrison et al. 2012:
%
T = 0.74 + 0.44;
L = 1.5;
durSwing = 0.44;
tStanceToSwing = T - durSwing;
tEMGdelay = 0.05; 
omegaEMG = .02; 
h = 0.5;
thEq = 0;

% %horse overground trot, Table 2 and Fig. 8 of Harrison et al. 2012:
%
% T = 0.29 + 0.39;
% L = 1;
% durSwing = 0.39;
% tStanceToSwing = T - durSwing;
% tEMGdelay = 0.05; 
% h = 0.5;


tStart = tStanceToSwing+tEMGdelay;
tEnd = 2*T+tEMGdelay;


hFig = figure;

subplot(3,3,1)
title('agonist (extensor) activation')
ylabel('a.u.')
colororder(winter(length(T)));

subplot(3,3,2)
title('agonist activation (t.n.)')
colororder(winter(length(T)));

subplot(3,3,3)
title('agonist activation (t. & a.n.)')
colororder(winter(length(T)));

subplot(3,3,4)
title('antag. (flexor) activation')
ylabel('a.u.')
colororder(winter(length(T)));

subplot(3,3,5)
title('antagonist activation (t.n.)')
colororder(winter(length(T)));

subplot(3,3,6)
title('antagonist activation (t. & a.n.)')
colororder(winter(length(T)));

subplot(3,3,7)
title('motion')
ylabel('rad')
colororder(winter(length(T)));
xlabel('time (s)')

subplot(3,3,8)
title('motion (t.n.)')
colororder(winter(length(T)));
xlabel('time (t.n.)')

subplot(3,3,9)
title('motion (t. & a.n.)')
colororder(winter(length(T)));
xlabel('time (t.n.)')

m = 3;
numT = length(T);


t = linspace(0,T,n+1)';

[theta,thDot,thDDot] = asymmetricalPolynomial(t,h,tStanceToSwing);


t(end) = [];
theta(end) = [];
thDot(end) = [];
thDDot(end) = [];

t3 = [];
for j=1:m
    t3 = [t3;(j-1)*T + t]; %#ok<AGROW>
end
theta3 = repmat(theta,m,1);
thDot3 = repmat(thDot,m,1);
thDDot3 = repmat(thDDot,m,1);

baseParams = [m0,c0,k0,s,a,thEq];
baseNoViscoelasticParams = [m0,0,0,s,a,thEq];
% baseNoViscoelasticParams = [0,c0,k0,s,0,thEq];

[agEMG3,antagEMG3,tEMG3] = approximateEMG(t3,theta3,tStanceToSwing,omegaEMG,T,L,baseParams);

figure(hFig);
hold on

bmTimeOfInterest = (t3 >= tStart) & (t3 <= tEnd);
p = sum(bmTimeOfInterest);
agEMG = agEMG3(bmTimeOfInterest);
antagEMG = antagEMG3(bmTimeOfInterest);
tEMG = tEMG3(bmTimeOfInterest);
theta = theta3(bmTimeOfInterest);
alpha = thDDot3(bmTimeOfInterest);

tshift = tStanceToSwing; %max(tEMG);

subplot(3,3,1)
area([T,T+tStanceToSwing]-tshift,max(max(agEMG),max(antagEMG))+[0,0],min(min(agEMG),min(antagEMG)),'facecolor',.9*[1,1,1],'edgealpha',0)
hold on
plot(tEMG-tshift,agEMG,'linewidth',lw)
drawnow
xlim([tStart,tEnd]-tshift)

subplot(3,3,2)
hold on
plot(tEMG/T,agEMG,'linewidth',lw)
% xlim([0,1])

subplot(3,3,3)
hold on
plot(tEMG/T,agEMG/max(agEMG),'linewidth',lw)
% xlim([0,1])

%     antagEMG = -min(0,EMG);

subplot(3,3,4)
area([T,T+tStanceToSwing]-tshift,max(max(agEMG),max(antagEMG))+[0,0],min(min(agEMG),min(antagEMG)),'facecolor',.9*[1,1,1],'edgealpha',0)
hold on
plot(tEMG-tshift,antagEMG,'linewidth',lw)
xlim([tStart,tEnd]-tshift)

subplot(3,3,5)
hold on
plot(tEMG/T,antagEMG,'linewidth',lw)
% xlim([0,1])

subplot(3,3,6)
hold on
plot(tEMG/T,antagEMG/max(antagEMG),'linewidth',lw)
% xlim([0,1])

subplot(3,3,7)
area([T,T+tStanceToSwing]-tshift,h+[0,0],-h,'facecolor',.9*[1,1,1],'edgealpha',0)
hold on
plot(tEMG-tshift,theta,'linewidth',lw)
%     if length(T) == 1
%         plot([0,tStanceToSwing]-tshift,1.5*min(theta)*[1,1],'k','linewidth',2)
%     end
xlim([tStart,tEnd]-tshift)
xlabel('time (s)')

subplot(3,3,8)
hold on
plot(tEMG/T,theta,'linewidth',lw)
% xlim([0,1])

subplot(3,3,9)
hold on
plot(tEMG/T,theta,'linewidth',lw)
% xlim([0,1])

drawnow

figure(hFig);
ax1 = subplot(3,3,1);
ax2 = subplot(3,3,4);

yLimMax = max(ax1.YLim(2),ax2.YLim(2));
ax1.YLim(2) = yLimMax;
ax2.YLim(2) = yLimMax;

[agEMGnovisc3,antagEMGnovisc3,tEMGnovisc3] = approximateEMG(t3,theta3,tStanceToSwing,omegaEMG,T,L,baseNoViscoelasticParams);

agEMGnovisc = agEMGnovisc3(bmTimeOfInterest);
antagEMGnovisc = antagEMGnovisc3(bmTimeOfInterest);
tEMGnovisc = tEMGnovisc3(bmTimeOfInterest);


h2 = figure;
sp1 = subplot(1,4,1);
area([T,T+tStanceToSwing]-tshift,.75*h+[0,0],-.75*h,'facecolor',.9*[1,1,1],'edgealpha',0)
hold on
set(sp1,'ColorOrderIndex',1);
plot(tEMG-tshift,theta-tEMGdelay,'linewidth',lw)
ylim([-.7*h,.7*h])
xlim([tStanceToSwing,2*T]-tshift)

sp2 = subplot(1,4,2);
hdd = max(alpha);
area([T,T+tStanceToSwing]-tshift,1.5*hdd+[0,0],-1.5*hdd,'facecolor',.9*[1,1,1],'edgealpha',0)
hold on
set(sp1,'ColorOrderIndex',1);
plot(tEMG-tshift,alpha-tEMGdelay,'linewidth',lw)
% ylim([-.7*h,.7*h])
xlim([tStanceToSwing,2*T]-tshift)

sp3 = subplot(1,4,3);
area([T,T+tStanceToSwing]-tshift,max(max(agEMG),max(antagEMG))+[0,0],min(min(agEMG),min(antagEMG)),'facecolor',.9*[1,1,1],'edgealpha',0)
hold on
set(sp3,'ColorOrderIndex',1);
plot(tEMG-tshift-tEMGdelay,agEMG,'linewidth',lw)
plot(tEMG-tshift-tEMGdelay,antagEMG,'linewidth',lw)
xlim([tStanceToSwing,2*T]-tshift)

sp4 = subplot(1,4,4);
area([T,T+tStanceToSwing]-tshift,max(max(agEMGnovisc),max(antagEMGnovisc))+[0,0],min(min(agEMGnovisc),min(antagEMGnovisc)),'facecolor',.9*[1,1,1],'edgealpha',0)
hold on
set(sp4,'ColorOrderIndex',1);
plot(tEMG-tshift-tEMGdelay,agEMGnovisc,'linewidth',lw)
plot(tEMG-tshift-tEMGdelay,antagEMGnovisc,'linewidth',lw)
xlim([tStanceToSwing,2*T]-tshift)

set(h2,'Position',[186.6000 265 1.0768e+03 205.6000])



