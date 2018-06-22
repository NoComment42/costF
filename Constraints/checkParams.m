clear;clc;

v = 0.5;
s = 2.0;
t = 1.0;
x = 1.4;
xp = abs(x-t);

% generalized exp
b = -2.0*v*s/log(2.0)
n = log(2.0)/v^b
g = exp(-n*xp^b);
eval = 1 - g

% inverse poly
b = -4*v*s
n = v^(-b)
g = 1/( n*xp^b + 1 );
eval = 1 - g
