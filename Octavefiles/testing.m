clear; clc;


xx = [ 1, 4, 6 ];
yx = [ 2, 6, 1 ];

a1 = (yx(2) - yx(1))/(xx(2) - xx(1))
b1 = yx(1) - a1*xx(1)

a2 = (yx(3) - yx(2))/(xx(3) - xx(2))
b2 = yx(2) - a2*xx(2)

x = [];
y1 = [];
y2 = [];
step = 0.5;
start = 1-step;
for i=1:6/step
  x(i) = i*step+start;
  y1(i) = piecewise1(x(i));
  y2(i) = piecewise2(x(i));
endfor
x
plot(x,y1,x,y2);