clc;
close all;
clear all;
i=1;
for D = 0:0.001:1
    Vin=24;
    C=10e-6;
    L=330e-6;
    R=20;
    RLesr=100e-3;
    RCesr=100e-3;
    Vd=0.7;
    Rdson=1e-3;
    RT=RLesr+Rdson;
    ko=1+RCesr/R;

    a11=-(1/L)*(  RT +(RCesr/ko)*(1-D)  );
    a12=-(1/L)*(  ( 1-RCesr/(ko*R) )*(1-D)  );
    a21= (1/C)*(  (1-D)/ko  );
    a22=-(1/C)*(  1/(ko*R)  );
    b1=  (Vin/L);

    A=[a11 a12;a21 a22];
    B=[b1;0];
    % x = [iL vC]^T

    xfix = -(A^(-1))*B;

    iLfix = xfix(1);
    vCfix = xfix(2);
    
    iC = iLfix*(1-D)/ko - vCfix/(ko*R);
    Dd(i) = D;
    effI(i) = (iC*RCesr +vCfix)^2/(R*Vin*iLfix);
    md(i)   = (iC*RCesr +vCfix)/Vin;
    i=i+1;
end
subplot(2,1,1)
plot(Dd,effI*100)
ylabel('$\eta=\frac{P_{out}}{R_{in}}$','Interpreter','latex','fontsize',24)
subplot(2,1,2)
plot(Dd,md)
ylabel('$M(D)=\frac{V_{out}}{V_{in}}$','Interpreter','latex','fontsize',24)
xlabel('$D$','Interpreter','latex','fontsize',24)