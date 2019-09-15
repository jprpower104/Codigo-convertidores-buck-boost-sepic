function evalBuck
clear;
clc;

tmax=0.002;
Ts=1e-7;
t=[0:Ts:tmax];

options=odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6]);
[T,Y] = ode45(@BuckSW,t,[0 0],options);

options=odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6]);
[T1,Y1] = ode45(@BuckPRO,t,[0 0],options);

subplot(2,1,1);
plot(T,Y(:,1),'b',T1,Y1(:,1),'r');
title('$Convertidor~Buck~Reductor~([0~1.0])$','Interpreter','latex','fontsize',15);
ylabel('$i_{L}(t)$','Interpreter','latex','fontsize',24);
subplot(2,1,2);
plot(T,Y(:,2),'b',T1,Y1(:,2),'r');
ylabel('$v_{C}(t)$','Interpreter','latex','fontsize',24);
xlabel('$t$','Interpreter','latex','fontsize',24);
end
function dy = BuckPRO(t,y)
dy = zeros(2,1);    % a column vector

T=1/100e3;
D=0.545;
Vin=22;
C=4.2e-6;
L=5e-6;
R=5;
K=(2*L)/(R*T);
Vc=2*Vin/(    1+sqrt( 1 + (4*K)/(D^2) )    );
D1=(Vin*D-Vc*D)/Vc;

A=[0 -(D+D1)/L;(D+D1)/C -1/(R*C)];
B=[Vin*D/L;0];

if (y(1,:)<0)
    y(1,:)=0;
end
dy=A*y+B;
end

function dy = BuckSW(t,y)
dy = zeros(2,1);    % a column vector

T=1/100e3;
D=0.545;
Vin=22;
C=4.2e-6;
L=4.32e-6;
R=5;

K=(2*L)/(R*T);
Vc=2*Vin/(   1+sqrt( 1+(4*K)/(D^2) )   );
D1=( Vin*D-Vc*D )/Vc;

h=(1/T)*mod(t,T);
u=compare(D,h);
u1=compare(D+D1,h);

if (u1==1) 
    A=[0 -1/L;1/C -1/(R*C)];
    B=[Vin*u/L;0];
else
    A=[0 0;0 -1/(R*C)];
    B=[Vin*u/L;0];
end

if (y(1,:)<0)
    y(1,:)=0;
end

dy=A*y+B;
end

function y2 = sat3(x2)
y2=1/2*(1+abs(x2)-abs(x2-1));
end
function a= compare(d,h2)
b=(tanh((h2-d)*1e12)+1)/2;
a=2-((tanh((b)*1e12)+1)/2)*2;
end