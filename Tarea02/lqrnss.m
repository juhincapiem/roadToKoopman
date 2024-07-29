%%%%%%%%%%%%%
%% The following is lqrnss.m
function [x,u,K]=lqrnss(As,Bs,Fs,Qs,Rs,xO,tspan)

global A E F Md tf W11 W12 W21 W22 n,

%
if nargin<7
error('Incorrect number of inputs specified')
return
end
%
%Convert Variables to normal symbology to prevent
% problems with global statement
%
A=As;
B=Bs;
F=Fs;
Q=Qs;
R=Rs;
plotflag=0; %set plotflag to 1 to avoid plotting of


[n,m]=size(A);
[nb,mb]=size(B);
[nq,mq] =size (Q) ;
[nr,mr]=size(R);
[nf ,mf] =size (F) ;

%
t0 = tspan(1);
tf = tspan(2);
tspan = [tf t0];
E=B*inv(R)*B';
Z=[A,-E;-Q,-A'];
[W, D] =eig (Z);

j=n;
[m1,index1]=sort(real(diag(D)));

for i=1:1:n
    m2(i)=m1(j);
    index2(i)=index1(j);
    index2(i+n)=index1(i+n);
    j=j-1;
end

Md=-diag(m2);

for i=1:2*n
    w2(:,i)=W(:,index2(i));
end
W=w2;

W11=zeros(n);
W12=zeros(n);
W21=zeros(n);
W22=zeros(n);

j=1 ;
for i=1:2*n:(2*n*n-2*n+1)
    W11(j:j+n-1)=W(i:i+n-1);
    W21(j:j+n-1)=W(i+n:i+2*n-1);
    W12(j:j+n-1)=W(2*n*n+i:2*n*n+i+n-1);
    W22(j:j+n-1)=W(2*n*n+i+n:2*n*n+i+2*n-1);
    j=j+n;
end

tl = 0.0;
tx = 0.0;
tu = 0.0;
x = 0.0;

[tx,x]=ode45('lqrnssf',fliplr(tspan),xO, ...
odeset('refine',2,'RelTol',le-4,'AbsTol',le-6));

j = 1;
us = 0.0; %Initialize computational variable
for i=1:1:mb
    for tua=t0:0.1:tf
        Tt=-inv(W22-F*W12)*(W21-F*W11);
        P=(W21+W22*expm(-Md*(tf-tua))*Tt* ...
        expm(-Md*(tf-tua)))*inv(W11+W12*expm(-Md*(tf-tua)) ...
        *Tt*expm(-Md*(tf-tua)));
        K=inv(R)*B'*P;
        xs=interp1(tx,x,tua);
        us1=real(-K*xs');
        us(j) =us1(i) ;
        tu(j)=tua;
        j=j+1;
    end
    u(:,i) =us' ;
    us=0;
    j=1;
end

%
%Provide final steady-state K
%
P=W21/Wll;
K=real(inv(R)*B'*P);
%
%Plotting Section, if desired
%
if plotflag~=l
%
%Plot diagonal Riccati coefficients using a
% flag variable to hold and change colors
%
fig=l ;
cflag=l;
j=1;
Ps=0.0;
%Figure number
%Variable used to change plot color
%Initialize P matrix plot variable
for i=1:1:n*n
    for tla= t0:0.1:tf
        Tt=-inv(W22-F*W12)*(W21-F*W11);
        P=real((W21+W22*expm(-Md*(tf-tla))*Tt*expm(-Md* ...
        (tf-tla)))*inv(Wll+W12*expm(-Md*(tf-tla))*Tt ...
        *expm(-Md*(tf-t1a))));
        Ps(j)=P(i);
        t1(j)=t1a;
        j=j+1;
    end
    if cflag==1;
        figure (fig)
        plot(tl,Ps, 'b')
        title('Plot of Riccati Coefficients')
        xlabel ('t')
        ylabel ("P Matrix")
        hold on;
        cflag=2;
    
    elseif cflag==2
        plot (t1, Ps, "m:")
        cflag=3;
    elseif cflag==3
        plot(t1,Ps,'g-.')
        cflag=4;
    elseif cflag==4
        plot(t1,Ps,'r--')
        cflag=1 ;
        fig=fig+1;
    end
Ps=0.0;
j=1;

end
if cflag==2 || cflag==3 || cflag==4
    hold
    fig=fig+1;
end
%
%Plot Optimized x
%
if n>2
    for i=1:3:(3*fix((n-3)/3)+1)
        Appendix C: MATLAB Files
        figure(fig);
        plot(tx,real(x(:,i)),'b',tx,real(x(:,i+1)),'m:',tx, ...
        real(x(:,i+2)),'g-.')

    end
    title('Plot of Optimized x')
    xlabel ("t")
    ylabel('x vectors')
    fig=fig+1;
end
if (n-3*fix(n/3))==1
    figure(fig);
    plot(tx,real(x(:,n)),'b')
elseif (n-3*fix(n/3))==2
    figure(fig);
    plot (tx,real(x(:,n-1) ) , "b" ,tx, real(x(:,n)),"m:");
end

title('Plot of Optimized x')
xlabel ('t')
ylabel('x vectors')
fig=fig+1;

if mb>2
    for i=1:3:(3*fix((mb-3)/3)+1)
    figure(fig);
    plot(tu,real(u(:,i)),'b',tu,real(u(:,i+1)),'m:', ...
    tu,real(u(:,i+2)),'g-.')
    title('Plot of Optimized u')
    xlabel ('t')
    ylabel('u vectors')
    fig=fig+1;
    end
end
if (mb-3*fix(mb/3))==1
    figure(fig);
    plot(tu,real(u(:,mb)),'b')
elseif (mb-3*fix(mb/3))==2
    figure(fig);
    plot(tu,real(u(:,mb-1)),'b',tu,real(u(:,mb)),'m:')
end
title('Plot of Optimized u')
xlabel (' t')
ylabel('u vectors')

end