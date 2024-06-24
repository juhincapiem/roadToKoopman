% Finite-Time Optimal Control
% Finite-Time LQR double integrator
% preambule
clear; close all; clc;
set(0,'DefaultLineLineWidth',2) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1]) %background white on plots
%dynamics
n=2; m=1; tf =30;
A = [0 1;0 0];
f = @(t,x) A*x;
B = [0;1];
g = @(t,x) B;
ft = @(t,x,u) f(t,x)+g(t,u)*u;
q =1;
Q_x = [q^2 0; 0 0];
Q_T = eye(n);
R =1;
%% control Desing
dtt= 0.1;
tSpan = 0:dtt:tf;  % Forward
%tSpan = fliplr(tSpan);  %Integrate backwards in time
options = odeset('RelTol',1e-13,'AbsTol',1e-300);
z0 = reshape(Q_T,n*n,1);
%for k=1:length(tSpan)-1
%for k=1:5
[tspan, P_i] = ode45(@(t,z)Mat2Vec(t,z,A,B,Q_x, R ,n),tSpan, z0, options);
% figures
plot(tspan, P_i(:,1),tspan, P_i(:,2),tspan, P_i(:,3),tspan, P_i(:,4))
%% simulation
%
%Runge-Kutta 4
k1 = @(t,x,u) (  ft(t,x,u) );
k2 = @(t,x,u) ( ft(t,x + k1(t,x,u)*dtt/2,u) );
k3 = @(t,x,u) ( ft(t,x + k2(t,x,u)*dtt/2,u) );
k4 = @(t,x,u) ( ft(t,x + k1(t,x,u)*dtt,u) );
f_ud = @(t,x,u) ( x + (dtt/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );
%
x0 =[-1;1];
x_ol = zeros(n,length(tSpan));
x_ol(:,1) = x0;
x_cl = zeros(n,length(tSpan));
x_cl(:,1) = x0;
for i = 1:length(tSpan)-1
    x_ol(:,i+1)=f_ud(0,x_ol(:,i),0);  %open loop
    u = -(R\B')*reshape(P_i(i,:),n,n)*x_cl(:,i);  % u=
    x_cl(:,i+1)=f_ud(0,x_cl(:,i),u);
end
figure
%plot(tSpan, x_ol(1,:), tSpan, x_ol(2,:));
plot(tSpan, x_cl(1,:), tSpan, x_cl(2,:));


% additional functions
%transform the matrix differential equation to vector ODE
function dz = Mat2Vec(~,z,A,B,Q,R,n)
P = reshape(z,n,n);
dP = (A'*P + P*A - P*B*(R\B'*P) + Q); %forward
%dP = -(A'*P + P*A - P*B*(R\B'*P) + Q); %backward
dz = reshape(dP,n*n,1);
end