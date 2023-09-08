
n = 500;
x = linspace(-pi/2,pi/2,2*n+1);
dx = x(2)-x(1);
lw = 1;

err = (x - sin(x))./x;

h1 = figure;
subplot(1,3,1)
area([-pi/6,pi/6],[max(err(abs(x)<=pi/3)),max(err(abs(x)<=pi/3))],'facecolor','k','facealpha',0.1,'edgecolor','none')
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(x,err,'linewidth',lw)
ylabel('Error in small angle assumption')
xlabel('Joint angle (degrees)')
xticks((-60:30:60)*pi/180);
xticklabels({'-60','-30','0','30','60'})
xlim([-pi/3,pi/3])
ylim([0,max(err(abs(x)<=pi/3))])
title('A.')
grid on

f = (sin(x)./x).^2;
% f = (sin(2*x)./(2*x)).^2;
f(n+1) = 1;

subplot(1,3,2)
area([-pi/6,pi/6],[max(f),max(f)],'facecolor','k','facealpha',0.1,'edgecolor','none')
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(x,f,'linewidth',lw)
ylabel('Apparent inertia scaling')
xlabel('Joint angle (degrees)')
xticks((-60:30:60)*pi/180);
xticklabels({'-60','-30','0','30','60'})
xlim([-pi/3,pi/3])
ylim([0,max(f)])
title('B.')
grid on

g = zeros(size(f));

for i=0:n
    g(n+1+i) = mean(f(n+1-i:n+1+i));
    g(n+1-i) = mean(f(n+1-i:n+1+i));
end

subplot(1,3,3)
area([0,pi/6],[max(g),max(g)],'facecolor','k','facealpha',0.1,'edgecolor','none')
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(x,g,'linewidth',lw)
ylabel('Average inertia error')
xlabel('Joint excursion (degrees)')
xticks((0:30:90)*pi/180);
xticklabels({'0','30','60','90'})
xlim([0,pi/2])
ylim([0,max(g)])
title('C.')
grid on

h1.Position(3) = 650;
h1.Position(4) = 200;


h2 = figure;
subplot(1,3,1)
area([-pi/6,pi/6],[max(err(abs(x)<=pi/3)),max(err(abs(x)<=pi/3))],'facecolor','k','facealpha',0.1,'edgecolor','none')
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(x,(x - sin(x))./x,'linewidth',lw)
ylabel('Error in small angle assumption')
xlabel('Joint angle (degrees)')
xticks((-60:30:60)*pi/180);
xticklabels({'-60','-30','0','30','60'})
xlim([-pi/3,pi/3])
ylim([0,max(err(abs(x)<=pi/3))])
title('A.')
grid on

f = (x./sin(x)).^2;
f(n+1) = 1;

subplot(1,3,2)
area([-pi/6,pi/6],[1.5,1.5],'facecolor','k','facealpha',0.1,'edgecolor','none')
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(x,f,'linewidth',lw)
ylabel('Apparent inertia scaling')
xlabel('Joint angle (degrees)')
xticks((-60:30:60)*pi/180);
xticklabels({'-60','-30','0','30','60'})
xlim([-pi/3,pi/3])
ylim([0,1.5])
title('B.')
grid on

g = zeros(size(f));

for i=0:n
    g(n+1+i) = mean(f(n+1-i:n+1+i));
    g(n+1-i) = mean(f(n+1-i:n+1+i));
end

subplot(1,3,3)
area([0,pi/6],[1.5,1.5],'facecolor','k','facealpha',0.1,'edgecolor','none')
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(x,g,'linewidth',lw)
ylabel('Average inertia error')
xlabel('Joint excursion (degrees)')
xticks((0:30:90)*pi/180);
xticklabels({'0','30','60','90'})
xlim([0,pi/2])
ylim([0,1.5])
title('C.')
grid on

h2.Position(3) = 650;
h2.Position(4) = 200;