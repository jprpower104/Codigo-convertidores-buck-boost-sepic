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
D=0.20834;
Vin=24;
C=81e-9;
L=528e-6;
R=10;

A=[0 -1/L;1/C -1/(R*C)];
B=[Vin*D/L;0];

dy=A*y+B;
end

function dy = BuckSW(t,y)
dy = zeros(2,1);    % a column vector

T=1/100e3;
D=0.20834;
Vin=24;
C=81e-9;
L=528e-6;
R=10;
h=(1/T)*mod(t,T);
u=compare(D,h);

A=[0 -1/L;1/C -1/(R*C)];
B=[Vin*u/L;0];

dy=A*y+B;
end

function y = sat3(x)
y=1/2*(1+abs(x)-abs(x-1));
end
function a= compare(d,h)
b=(tanh((h-d)*1e12)+1)/2;
a=2-((tanh((b)*1e12)+1)/2)*2;
end