
%human hip, Stein, Zehr 1996
%human tibia, Stein, Zehr 1996
%human medial gastrocnemius, Maisetti, Hug, et al. 2012, approx. from Fig 2A.
%Hajian and Howe 1997, Fig. 6, with r = 9.76 cm (human finger length)
%stick insect FTi, von Twickel, Guschlbauer et al. 2019, Fig. 1c.
%Zakotnik, Matheson, and Duerr 2006, Fig. 8

L = [1,0.50,.22,.0976,4e-2,.0228]; %m
k = [100/(pi/2),(18.5-13.6),10/(pi/6),50*(.0976)^2, 4.5e-6/(2*pi/3) 90e-6/(pi/3)]; %Nm/rad

L(end-1) = [];
k(end-1) = [];

figure
plot(L,k,'o-')
ax = gca;
ax.YScale = 'log';
ax.XScale = 'log';
grid on
xlabel('L (m)')
ylabel('k_{r,elas} (Nm/rad)')
xlim([1e-2,1])

A = [log10(L)',1+zeros(size(L'))];
b = log10(k)';

params = linsolve(A,b);

scaleExp = params(1)

hold on
plot(L,10.^(params(2))*L.^params(1))