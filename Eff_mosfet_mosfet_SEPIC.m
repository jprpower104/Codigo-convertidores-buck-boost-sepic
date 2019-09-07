clc;
close all;
clear all;
i=1;
for D = 0:0.001:1
    Vin=20;

RL1=1e-3;
RC1=1e-3;
RL2=1e-3;
RC2=1e-3;

Ron1= 10e-3;
Ron2= 10e-3;

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
    % x = [iL vC]^T

    xfix = -(A^(-1))*B;

    iL1fix = xfix(1);
    iL2fix = xfix(2);
    vC1fix = xfix(3);
    vC2fix = xfix(4);
    
    iC2 = (iL1fix-iL2fix)*(1-D)/k0 - vC2fix/(k0*R0);
    Dd(i) = D;
    effI(i) = (iC2*RC2 +vC2fix)^2/(R0*Vin*iL1fix);
    md(i)   = (iC2*RC2 +vC2fix)/Vin;
    i=i+1;
end
subplot(2,1,1)
plot(Dd,effI*100)
ylabel('$\eta=\frac{P_{out}}{R_{in}}$','Interpreter','latex','fontsize',24)
subplot(2,1,2)
plot(Dd,md)
ylabel('$M(D)=\frac{V_{out}}{V_{in}}$','Interpreter','latex','fontsize',24)
xlabel('$D$','Interpreter','latex','fontsize',24)