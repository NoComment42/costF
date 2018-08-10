clear;clc;

data = load("../PRINTFILES/cp_constraint.dat");
data2 = load("../PRINTFILES/cf_constraint.dat");
data3 = load("../PRINTFILES/sep_constraint.dat");
data4 = load("../PRINTFILES/attach_constraint.dat");

subplot(2,2,1);
plot(data(:,1),data(:,2));

subplot(2,2,2);
plot(data2(:,1),data2(:,2));

subplot(2,2,3);
plot(data3(:,1),data3(:,2));

subplot(2,2,4);
plot(data4(:,1),data4(:,2));