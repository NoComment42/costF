clear; clc;clf;

f = load('functions.dat');
x = f(:,1);
f1 = f(:,2);
f2 = f(:,3);
f3 = f(:,4);

plot(x,f1,x,f2,x,f3);