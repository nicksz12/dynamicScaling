
%Nicholas Szczecinski
%West Virginia University
%13 October 2020
%Updated 8 September 2023

    
clearvars
close all

tPowerMin = -3;
tPowerMax = 2;
lPowerMin = -5.5;
lPowerMax = 1;

m0 = 12;
c0 = 1.31e3; %see zetaAcrossScales.m
k0 = 12e3; 
s = sqrt(1e-3);

baseParams = [m0,c0,k0,s];

Tn = 2*pi*sqrt(m0/k0);

%List of labels, characteristic times, and characteristic lengths.
%Name, period of locomotion (or other behavior), limb length, t-axis label placement, L-axis label placement.
            
behaviorCellSwing = {'fruit fly',            [.95, 1.05]*60e-3,    0.5e-3, 1.5,    1;... %Wosnitza et al. 2013
                'mouse',                [.95, 1.05]*.1,          2e-2,   1.5,    1;... %Herbin et al. 2007
                'rat',                  [.95, 1.05]*.2,           5e-2,   1.5,    1;... %Hruska et al. 1979
                'cat',                  [.95, 1.05]*250e-3,   25e-2,  1.5,    1;... %Grillner 1975
                'human',                [.95, 1.05]*500e-3,        1,      1.5,    1;... %Grillner et al. 1979
                'stick insect',         [.95, 1.05]*.6,          4e-2,   1.5,    1;... %Cruse and Bartling 1995
                'American cockroach',   [.95, 1.05]*.0435,     8e-3,   1.5,    1;... %Delcomyn 1971, The Locomotion of the Cockroach Periplaneta Americana
                'horse',                [.95, 1.05]*.44,          1.5,    1.5,    1;... %Hildebrand 1959, Hooper et al. 2009
                'tardigrade',           [.96, 1.05]*.25,        10e-6,  1.5,    1}; %Nirody, Duran, Johnston, and Cohen 2021
            
behaviorCellStance = {'fruit fly',            [60e-3, 250e-3],    0.5e-3, 1.5,    1;... %Wosnitza et al. 2013
                'mouse',                [.1, .33],          2e-2,   1.5,    1;... %Herbin et al. 2007
                'rat',                  [.2, .6],           5e-2,   1.5,    1;... %Hruska et al. 1979
                'cat',                  [250e-3, 700e-3],   25e-2,  1.5,    1;... %Grillner 1975
                'human',                [500e-3, 1],        1,      1.5,    1;... %Grillner et al. 1979
                'stick insect',         [.6, 1.8],          4e-2,   1.5,    1;... %Cruse and Bartling 1995
                'American cockroach',   [.0435, 0.667],     8e-3,   1.5,    1;... %Delcomyn 1971, The Locomotion of the Cockroach Periplaneta Americana
                'horse',                [.44, .9],          1.5,    1.5,    1;... %Hildebrand 1959, Hooper et al. 2009
                'tardigrade',           [0.5, 2.5],         10e-6,  1.5,    1}; %Nirody, Duran, Johnston, and Cohen 2021

            
nSampsPlot = 200;
nSampsInvDyn = 5;

%Plot how the forces that contribute to the steady-state response scale
%with L and T. This plot is a lot to take in, but it gives a useful
%summary.
hSteadyStateResponse = scalingInverseDynamicsAcrossScales(tPowerMin,tPowerMax,lPowerMin,lPowerMax,baseParams,nSampsInvDyn);
hSteadyStateResponseStance = scalingInverseDynamicsAcrossScalesStance(tPowerMin,tPowerMax,lPowerMin,lPowerMax,baseParams,nSampsInvDyn);

%Plot how the phase lag of the limb displacement, phi, scales with L and T.
%This is a proxy for where the actuator energy goes during steady-state
%motion.
%Plot how the damping ratio zeta scales with L and T. This is a proxy for 
%where unwanted energy from perturbations goes.
%Plot the size and cycle duration of various animal behaviors on the L and
%T plot. This shows that animal behaviors span several dynamic regimes, and
%some animals operate within several.
[fig2swing,fig3swing,fig4swing,figSswing] = phiAndX(tPowerMin,tPowerMax,lPowerMin,lPowerMax,behaviorCellSwing,baseParams,nSampsPlot);
[fig2stance,fig3stance,fig4stance,figSstance] = phiAndXstance(tPowerMin,tPowerMax,lPowerMin,lPowerMax,behaviorCellStance,baseParams,nSampsPlot);

%Plot simulated system responses. The simulateJointResponse
%function is useful for this.
%QS, overdamped
L = 1e-2;
T = 1;
[~,~,~,~,~,~,~,~,h] = simulateJointResponse(k0,c0,m0,s,L,1,T,3*T,[],[]);
figure(h);
sgtitle('Quasi-static, overdamped')
%QS, underdamped
L = 1.7;
T = 10;
[~,~,~,~,~,~,~,~,h] = simulateJointResponse(k0,c0,m0,s,L,1,T,3*T,[],[]);
figure(h);
sgtitle('Quasi-static, underdamped')
%Inertial, underdamped
L = 1.7;
T = 0.5;
[~,~,~,~,~,~,~,~,h] = simulateJointResponse(k0,c0,m0,s,L,1,T,3*T,[],[]);
figure(h);
sgtitle('Inertial, underdamped')
%viscous, overdamped
L = 1e-2;
T = 1e-2;
[~,~,~,~,~,~,~,~,h] = simulateJointResponse(k0,c0,m0,s,L,1,T,3*T,[],[]);
figure(h);
sgtitle('Viscous, overdamped')
%Inertial, overdamped
L = .03;
T = 10^-2.8;
[~,~,~,~,~,~,~,~,h] = simulateJointResponse(k0,c0,m0,s,L,1,T,3*T,[],[]);
figure(h);
sgtitle('Inertial, overdamped')
