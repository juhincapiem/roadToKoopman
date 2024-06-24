% Finite-Time Optimal Control 
% We simulate a scalar analytical example
% preambule
clear; close all; clc;
set(0,'DefaultLineLineWidth',2) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1]) %background white on plots
%dynmics
a = 3; b=1; c=1; t0=0; tf=6;
f = @(t,x) a*x;
g = @(t,x) b;
ft = @(t,x,u) f(t,x)+g(t,u)*u;

%ctrl design
u_opt = @(t,x0) -( 2*a*b*c*exp( a*(2*tf-t0-t) )*x0 ) / ( 2*a-b^2*c*(1-exp(2*a*(tf-t0))) );
x_state = @(t,x0) x0*exp( a*(t-t0) ) + ((b^2*c*exp(a*(tf-t0))*x0)/(2*a-b^2*c*(1-exp(2*a*(tf-t0)))))...
                   *(exp( a*(tf-t) ) -exp( a*(t+tf-2*t0) ));
%% simulation
dtt= 0.1;
tspan = 0:dtt:tf;
%tspan = [t0 tf];
x0=5;
%options = odeset('RelTol',1e-13,'AbsTol',1e-300);
options = odeset('RelTol',1e-13,'AbsTol',1e-300);
[tspan, x_OpLoop] = ode45(@(t,x)ft(t, x, 0), tspan, x0, options);
[tspan1, x_opt] = ode23s(@(t,x)ft(t, x, u_opt(t,x0)), tspan, x0, options);
% Figures
%plot(tspan, x_OpLoop,'--r'); 
plot(tspan1, x_opt,'--r'); hold on
plot(tspan1, x_state(tspan1,x0),'k'); 
grid on
%title('LQR')
xlabel('$t$ ','interpreter','latex');
ylabel('$x$','interpreter','latex');
set(gca,'fontsize',20)
LEG = legend('$x$ ');
set(LEG,'interpreter','latex')