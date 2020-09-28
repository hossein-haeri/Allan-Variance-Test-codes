clear all
clc
close all

x = linspace(0,4*pi,10);
y = sin(x);
p = polyfit(x,y,3);

x1 = linspace(0,4*pi);
y1 = polyval(p,x1);
figure
plot(x,y,'o')
hold on
plot(x1,y1)
hold off