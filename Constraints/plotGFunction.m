clear;clc;

data = load("c1.dat");
data2 = load("c2.dat");

figure(1);
plot(data(:,1),data(:,2));
figure(2);

plot(data2(:,1),data2(:,2));