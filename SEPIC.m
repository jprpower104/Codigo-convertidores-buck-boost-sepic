function SEPIC
clear;
clc;

tmax=0.1;
Ts=1e-7;
t=[0:Ts:tmax];
error=1e-12;
options=odeset('RelTol',error);%paso los vectores de 2 a 4
[T,Y] = ode23s(@SepicSW,t,[0 0 0 0],options);                 %paso los vectores de 2 a 4

options=odeset('RelTol',error);%paso los vectores de 2 a 4
[T1,Y1] = ode45(@SepicPRO,t,[0 0 0 0],options);              %paso los vectores de 2 a 4

subplot(4,1,1);
plot(T,Y(:,1),'b',T1,Y1(:,1),'r');
title('$Convertidor~Sepic~([0~1.0])$','Interpreter','latex','fontsize',15);
ylabel('$i_{L_1}(t)$','Interpreter','latex','fontsize',24);
subplot(4,1,2);
plot(T,Y(:,2),'b',T1,Y1(:,2),'r');
ylabel('$i_{L_2}(t)$','Interpreter','latex','fontsize',24);
subplot(4,1,3);
plot(T,Y(:,3),'b',T1,Y1(:,3),'r');
ylabel('$v_{C_1}(t)$','Interpreter','latex','fontsize',24);
subplot(4,1,4);
plot(T,Y(:,4),'b',T1,Y1(:,4),'r');
ylabel('$v_{C_2}(t)$','Interpreter','latex','fontsize',24);
xlabel('$t$','Interpreter','latex','fontsize',24);
end
function dy = SepicPRO(t,y)
dy = zeros(4,1);    % a column vector

D=0.8;
Vin=12;
C1=32e-6;
C2=8e-6;
L1=50e-6;
L2=200e-6;
R=10;

A=[0        0    -(1-D)/L1   -(1-D)/L1 ;
   0        0       -D/L2     (1-D)/L2 ;   
 (1-D)/C1   D/C1      0          0     ;      
 (1-D)/C2  -(1-D)/C2  0      -1/(R*C2)];

B=[Vin/L1;0;0;0];

dy=A*y+B;
end

function dy = SepicSW(t,y)
dy = zeros(4,1);    % a column vector

T=1/100e3;

D=0.8;
Vin=12;
C1=32e-6;
C2=8e-6;
L1=50e-6;
L2=200e-6;
R=10;

h=(1/T)*mod(t,T);
u=compare(D,h);
 
A=[0        0    -(1-u)/L1  -(1-u)/L1 ;
   0        0      -u/L2     (1-u)/L2 ;   
 (1-u)/C1   u/C1      0         0     ;      
 (1-u)/C2 -(1-u)/C2   0     -1/(R*C2)];

B=[Vin/L1;0;0;0];

dy=A*y+B;
end


%Verificar de aca para abajo
function y = sat3(x)
y=1/2*(1+abs(x)-abs(x-1));
end
function a= compare(d,h)
b=(tanh((h-d)*1e12)+1)/2;
a=2-((tanh((b)*1e12)+1)/2)*2;
end