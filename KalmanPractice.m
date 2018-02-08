close all;
clear all;

%scalar Kalman filter

%model:
%x(t+1) = Ax(t)+Bu(t)+v(t)
%y(t) = Cx(t)+w(t)
std_dev = 2;
%time
t = 0:0.1:20;
%covariance matrix of v:
V = std_dev*0.5;
%covariance matrix of w:
W = std_dev*std_dev;

%initial x:
x0 = 0;
P0 = 1; %variance
xpriori = x0;
Ppriori = P0;
xposteriori = x0;
Pposteriori = P0;
%matrices
A = 1;
B = 0;
C = 1;
K = 0;
S = 0;

y = zeros(size(t));
y_pred = zeros(size(t));
y_pred(1) = x0;
u = zeros(size(t));

for i=1:size(t,2)
    if i<100
        y(i) = std_dev*randn(1);
    else
        y(i) = std_dev*randn(1)+20;
    end
    
    if i>1
        %time update
        xpriori = A*xposteriori + B*u(i);
        Ppriori = A*Pposteriori*transpose(A) + V;
        %measurmenet update
        e = y(i) - C*xpriori;
        S = C*Ppriori*transpose(C) + W;
        K = Ppriori*transpose(C)*(S^(-1));
        xposteriori = xpriori + K*e;
        Pposteriori = Ppriori - K*S*transpose(K);
    
        y_pred(i) = xposteriori;
    end
end

plot(t,y, 'Color', 'red');
hold on;
plot(t,y_pred, 'Color', 'blue');
title('Kalman filter example and practice');
grid on;
legend('Measured value', 'Estimated value');
