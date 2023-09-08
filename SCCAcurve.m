function [y, v, a, j] = SCCAcurve(t,b,c,d,Ca)
    
    if b + c + d ~= 1
        warning('b, c, and d must sum to 1; setting c = 1 - b - d.')
        c = 1 - b - d;
    end

    x = t/t(end);

    y1 = Ca*(b/pi*x-(b/pi)^2*sin(pi/b*x));
    v1 = Ca*(b/pi-b/pi*cos(pi/b*x));
    a1 = Ca*sin(pi/b*x);
    j1 = Ca*pi/b*cos(pi/b*x);
    
    y2 = Ca*(x.^2/2+b*(1/pi-1/2)*x+b^2*(1/8-pi^-2));
    v2 = Ca*(x+b*(1/pi-1/2));
    a2 = Ca + zeros(size(x));
    j2 = zeros(size(x));
    
    y3 = Ca*((b/pi+c/2)*x+(d/pi)^2+b^2*(1/8-pi^-2)-(1-d)^2/8-(d/pi)^2*cos(pi/d*(x-(1-d)/2)));
    v3 = Ca*(b/pi+c/2+d/pi*sin(pi/d*(x-(1-d)/2)));
    a3 = Ca*cos(pi/d*(x-(1-d)/2));
    j3 = -Ca*pi/d*sin(pi/d*(x-(1-d)/2));
    
    y4 = Ca*(-x.^2/2+(b/pi+1-b/2)*x+(2*d^2-b^2)*(pi^-2-1/8)-1/4);
    v4 = Ca*(-x+b/pi+1-b/2);
    a4 = -Ca + zeros(size(x));
    j4 = zeros(size(x));
    
    y5 = Ca*(b/pi*x+2*(d^2-b^2)*pi^-2+((1-b)^2-d^2)/4-(b/pi)^2*sin(pi/b*(x-1)));
    v5 = Ca*(b/pi-b/pi*cos(pi/b*(x-1)));
    a5 = Ca*sin(pi/b*(x-1));
    j5 = Ca*pi/b*cos(pi/b*(x-1));
    
    %zone1 = (x >= 0 & x < b/2);
    zone2 = (x >= b/2 & x < (1-d)/2);
    zone3 = (x >= (1-d)/2 & x < (1+d)/2);
    zone4 = (x >= (1+d)/2 & x < 1-b/2);
    zone5 = (x >= 1-b/2);
        
    y = y1;
    y(zone2) = y2(zone2);
    y(zone3) = y3(zone3);
    y(zone4) = y4(zone4);
    y(zone5) = y5(zone5);
    
    v = v1;
    v(zone2) = v2(zone2);
    v(zone3) = v3(zone3);
    v(zone4) = v4(zone4);
    v(zone5) = v5(zone5);
    
    a = a1;
    a(zone2) = a2(zone2);
    a(zone3) = a3(zone3);
    a(zone4) = a4(zone4);
    a(zone5) = a5(zone5);
    
    j = j1;
    j(zone2) = j2(zone2);
    j(zone3) = j3(zone3);
    j(zone4) = j4(zone4);
    j(zone5) = j5(zone5);
    
    figure
    plot(t,y)
    hold on
    plot(t,v)
    plot(t,a)
    plot(t,j)
    
end