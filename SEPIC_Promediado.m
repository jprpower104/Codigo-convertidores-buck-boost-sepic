function SEPIC_Promediado
clear;
clc;

tmax=0.1;
Ts=1e-7;
t=[0:Ts:tmax];
error=1e-12;

options=odeset('RelTol',error);%paso los vectores de 2 a 4
[T,Y] = ode45(@SepicPRO,t,[0 0 0 0],options);              %paso los vectores de 2 a 4

subplot(4,1,1);
plot(T,Y(:,1),'b');
title('$Convertidor~Sepic~([0~1.0])$','Interpreter','latex','fontsize',15);
ylabel('$i_{L_1}(t)$','Interpreter','latex','fontsize',24);
subplot(4,1,2);
plot(T,Y(:,2),'b');
ylabel('$i_{L_2}(t)$','Interpreter','latex','fontsize',24);
subplot(4,1,3);
plot(T,Y(:,3),'b');
ylabel('$v_{C_1}(t)$','Interpreter','latex','fontsize',24);
subplot(4,1,4);
plot(T,Y(:,4),'b');
ylabel('$v_{C_2}(t)$','Interpreter','latex','fontsize',24);
xlabel('$t$','Interpreter','latex','fontsize',24);
end
function dy = SepicPRO(t,y)
dy = zeros(4,1);    % a column vector

D=0.8;
Vin=12;
C1=32e-6;
C2=32e-6;
L1=200e-6;
L2=200e-6;
R=10;

A=[0        0    -(1-D)/L1   -(1-D)/L1 ;
   0        0       -(D)/L2     (1-D)/L2 ;   
 (1-D)/C1   D/C1      0          0     ;      
 (1-D)/C2  -(1-D)/C2  0      -1/(R*C2)];

B=[Vin/L1;0;0;0];
 
dy=A*y+B;
%dy(1)=(-(1-D)/L1)*y(3)+(-(1-D)/L1)*y(4)+Vin/L1;
%dy(2)=(-(D)/L2)*y(3)+((1-D)/L2)*y(4);
%dy(3)=((1-D)/C1)*y(1)+(D/C1)*y(2);
%dy(4)=((1-D)/C2)*y(1)-((1-D)/C2)*y(2)-(1/(R*C2))*y(4);
end