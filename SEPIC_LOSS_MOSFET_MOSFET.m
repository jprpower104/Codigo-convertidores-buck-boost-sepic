function evalSEPIC_LOSS
clear;
clc;

tmax=0.01;
Ts=1e-6;
t=[0:Ts:tmax];
error =1e-6;
options=odeset('RelTol',error,'AbsTol',[error error error error]);
[T,Y] = ode45(@SEPIC_LOSS_SW,t,[0 0 0 0],options);

options=odeset('RelTol',error,'AbsTol',[error error error error]);
[T1,Y1] = ode45(@SEPIC_LOSS_PRO,t,[0 0 0 0],options);

subplot(4,1,1);
plot(T,Y(:,1),'b',T1,Y1(:,1),'r');
title('$Convertidor~SEPIC~([0~1.0])$','Interpreter','latex','fontsize',15);
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
function dy = SEPIC_LOSS_PRO(t,y)
dy = zeros(4,1);    % a column vector

D=0.7;
Vin=20;

RL1=1e-3;
RC1=2.3e-3;
RL2=2.5e-3;
RC2=3e-3;

Ron1= 1e-3;
Ron2= 2.5e-3;

L1=330e-6;
C1=2.2e-6;
L2=330e-6;
C2=22e-6;
R0=10;

k0=1+RC2/R0;

A=[-(1/L1)*( RL1+Ron1*D+(RC2/k0 + Ron2 + RC1)*(1-D) )   -(1/L1)*( Ron1*D + (Ron2+ RC2/k0)*(1-D) )             -(1/L1)*(1-D)     -(1/L1)*( 1-RC2/(k0*R0) )*(1-D);...
   -(1/L2)*( Ron1*D +(Ron2 + RC2/k0)*(1-D) )            -(1/L2)*( RL2+(RC1+Ron1)*D + (Ron2 + RC2/k0)*(1-D) )  -D/L2             (1/L2)*( 1-RC2/(k0*R0) )*(1-D);...
   (1-D)/C1                                             D/C1                                                  0                 0                            ;...
   (1-D)/(k0*C2)                                        -(1-D)/(k0*C2)                                        0                 -(1/C2)*( 1/(k0*R0) )];

B=[Vin/L1;...
   0;...
   0;...
   0];

dy=A*y+B;
end

function dy = SEPIC_LOSS_SW(t,y)
dy = zeros(2,1);    % a column vector

T=1/100e3;

D=0.7;
Vin=20;

RL1=1e-3;
RC1=2.3e-3;
RL2=2.5e-3;
RC2=3e-3;

Ron1= 1e-3;
Ron2= 2.5e-3;

L1=330e-6;
C1=2.2e-6;
L2=330e-6;
C2=22e-6;
R0=10;

k0=1+RC2/R0;


h=(1/T)*mod(t,T);
u=compare(D,h);

A=[-(1/L1)*( RL1 + Ron1*u + (RC2/k0 + Ron2 + RC1)*(1-u) )   -(1/L1)*( Ron1*u +( Ron2 + RC2/k0 )*(1-u) )             -(1/L1)*(1-u)   -(1/L1)*( 1-RC2/(k0*R0) )*(1-u);...
   -(1/L2)*( Ron1*u - (Ron2 + RC2/k0)*(1-u) )               -(1/L2)*( RL2 + (RC1+Ron1)*u + (Ron2 + RC2/k0)*(1-u) )  -u/L2           (1/L2)*( 1-RC2/(k0*R0) )*(1-u);...
    (1-u)/C1                                                u/C1                                                    0               0                            ;...
    (1-u)/(k0*C2)                                           -(1-u)/(k0*C2)                                          0               -(1/C2)*(1/(k0*R0))];

B=[Vin/L1;...
   0;...
   0;...
   0];

dy=A*y+B;
end

function y = sat3(x)
y=1/2*(1+abs(x)-abs(x-1));
end
function a= compare(d,h)
b=(tanh((h-d)*1e12)+1)/2;
a=2-((tanh((b)*1e12)+1)/2)*2;
end