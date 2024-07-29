function [t,Y,Y_m] = BDF_Diff_Riccati_Sym(A,B,C,Y0,t0,tf)
%  [t,Y,Y_m] = BDF_Diff_Riccati_Sym(A,B,C,Y0,t0,tf)
% The resolution of the symmetric continuous-time differential matrix
% Riccati equations ( backward differentiation formula method)
%
%                           dY/dt = A'Y + YA - YBB'Y + C'C (*)
%                           Y(t0) = Y0
%
% Input:        A          : size matrix  (m,m).
%               B          : size matrix  (m,s).
%               C          : size matrix  (p,m).
%               Y0         : symetric size matrix  (m,m).
%               tf , t0    : temps.
% Output:    
%             Y_m          : solution to the inslant Tf of (*)  
%             Y            : solution to (*) 
%             t
%    Author                :  LAKHLIFA SADEK.
%    Last modification     :  08/01/2019
% E-mail: lakhlifasdek@gmail.com; sadek.l@ucd.ac.ma
% ORCID : https://orcid.org/0000-0001-9780-2592
% 
if norm(Y0-Y0','fro')>1e-12 
    error('The coefficient matrix must be symmetric')
end
h=0.001;
N=(tf-t0)/h; 
[ms]=size(A,1);
Z0=Y0;Z1=Y0;
TAA=(2/3)*h*A-(1/2)*eye(ms); 
BB=sqrt((2/3)*h)*(B);
  alpha0=4/3;alpha1=-1/3;
t(1)=t0;
Y=[Y0(:)];
for k=1:N
    t(k+1)=t0+(k)*h;
    CC=alpha0*Z0+alpha1*Z1+h*(2/3)*(C'*C);
    [Yk,l,g]=care(TAA,BB,CC);
   Z1=Z0;
   Z0=Yk;  
   Y=[Y Yk(:)];
end
Y_m=Yk;
end