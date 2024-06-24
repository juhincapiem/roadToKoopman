% Código para resolver la ecuación de Burgers 1D
%           du/dt + u*du/dx = v*d^2u/dx^2
% Sujeto a las siguientes condiciones de frontera e iniciales
%       u(0,t) = u(L,t) = 0
%       u(x,0) = f(x)
% Código desarrollado por Juan Sebastián Hincapié para el proyecto
% de operador lineal Koopman en borde finito

% Propiedades de sistema
L = 2; %[m]
v = 0.08; %[m2/s]
tf = 1.0; %[s]

% Discretización del dominio temporal y espacial
Nt = 150;
Nx = 100;
x = linspace(0,L, Nx);
time = linspace(0,tf, Nt);
dx = x(2)-x(1);
dt = time(2)-time(1);

% Constantes del problemas
alfa = (dt*v)/(2*dx^2);
beta = dt/(2*dx);

% Matrices y vectores
A = zeros(Nx,Nx);
b = zeros(Nx,1);
campoVel = zeros(Nx,Nt);
uf = ones(Nx,1);
up = ones(Nx,1);
I = (1:1:Nx);
fig = 0;

% Llenamos la condición inicial del problema
for i=I(2:end-1)
    up(i) = U0(x(i));
end
campoVel(:,1) = up(:,1);

%Alertamos por el número de Fourir y por el número de Courant
sprintf("El máximo número CFL es de %0.5f",max(up)*dt/dx)
sprintf("El número de Fourier es de %0.5f",alfa*2)

% Verificamos las condiciones iniciales para la velocidad
fig = fig + 1;
figure(fig)
plot(x,up,"--r", LineWidth=1.2)
title("Condiciones iniciales")
ylabel("u [m/s]")
xlabel("x [m]")
grid on

% Resolvemos el problema para todo el tiempo
for t = 1:length(time)
    A(:,:) = matrixA(A,I(2:end-1),alfa,beta,up);
    b(:,1) = vectorB(b,I(2:end-1),alfa,beta,up);
    w(:,1) = control(x,t,I(2:end-1));
    campoVel(:,t) = up(:,1);
    uf = A\(b+w);
    up = uf;
end

% Verificamos las condiciones iniciales para la velocidad
fig = fig + 1;
figure(fig)
plot(x,up,"--r", LineWidth=1.2)
title(sprintf("Estado del sistema dinámico en t = %0.2f",tf))
ylabel("u [m/s]")
xlabel("x [m]")
grid on

%Creamos una animacion
animacion(campoVel,x,time)
