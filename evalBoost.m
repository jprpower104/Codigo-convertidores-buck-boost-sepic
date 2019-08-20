function evalBoost
clear;
clc;

tmax=4e-2;
Ts=1e-7;
t=[0:Ts:tmax];
error=1e-9;

options=odeset('RelTol',error,'AbsTol',[error error]);
[T,Y] = ode113(@BoostSW,t,[0 0],options);

options=odeset('RelTol',error,'AbsTol',[error error]);
[T1,Y1] = ode45(@BoostPRO,t,[0 0],options);

subplot(2,1,1);
plot(T,Y(:,1),'b',T1,Y1(:,1),'r');
title('$Convertidor~Boost~Elevador([0~1.0])$','Interpreter','latex','fontsize',15);
ylabel('$i_{L}(t)$','Interpreter','latex','fontsize',24);
subplot(2,1,2);
plot(T,Y(:,2),'b',T1,Y1(:,2),'r');
ylabel('$v_{C}(t)$','Interpreter','latex','fontsize',24);
xlabel('$t$','Interpreter','latex','fontsize',24);
end
function dy = BoostPRO(t,y)
dy = zeros(2,1);    % a column vector

d=0.6;%0.5
Vin=20;%24V
C=22e-6;%10u
L=330e-6;
R=100;

A=[0 -(1-d)/L;(1-d)/C -1/(R*C)];
B=[Vin/L;0];

dy=A*y+B;
% dy(1)=Vin/L-(1-d)*y(2)/L;
% dy(2)=(1-d)*y(1)/C-y(2)/(R*C);
end

function dy = BoostSW(t,y)
dy = zeros(2,1);    % a column vector

T=1/100e3;
d=0.6;%0.5
Vin=20;%24V
C=22e-6;%10u
L=330e-6;
R=100;

h=(1/T)*mod(t,T);
u=compare(d,h);

A=[0 -(1-u)/L;(1-u)/C -1/(R*C)];
B=[Vin/L;0];

dy=A*y+B;
% dy(1)=Vin/L-(1-u)*y(2)/L;
% dy(2)=(1-u)*y(1)/C-y(2)/(R*C);
end

function y = sat3(x)
y=1/2*(1+abs(x)-abs(x-1));
end
function a= compare(d,h)
b=(tanh((h-d)*1e12)+1)/2;
a=2-((tanh((b)*1e12)+1)/2)*2;
end