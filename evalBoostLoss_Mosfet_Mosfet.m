function evalBoostLoss
 cd
tmax=0.005;
Ts=1e-7;
t=[0:Ts:tmax];

options=odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6]);
[T,Y] = ode45(@BoostSWLoss,t,[0 0],options);

options=odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6]);
[T1,Y1] = ode45(@BoostPROLoss,t,[0 0],options);

subplot(2,1,1);
plot(T,Y(:,1),'b',T1,Y1(:,1),'r');hold on
title('$Convertidor~Boost~Elevador([0~1.0])$','Interpreter','latex','fontsize',15);
ylabel('$i_{L}(t)$','Interpreter','latex','fontsize',24);
subplot(2,1,2);
plot(T,Y(:,2),'b',T1,Y1(:,2),'r');hold on
ylabel('$v_{C}(t)$','Interpreter','latex','fontsize',24);
xlabel('$t$','Interpreter','latex','fontsize',24);
end

function dy = BoostPROLoss(t,y)
dy = zeros(2,1);    % a column vector

D=0.5;
Vin=24;
C=10e-6;
L=330e-6;
R=20;
RLesr=1e-3;
RCesr=1e-3;
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

dy=A*y+B;
end

function dy = BoostSWLoss(t,y)
dy = zeros(2,1);    % a column vector

T=1/100e3;
D=0.5;
Vin=24;
C=10e-6;
L=330e-6;
R=20;
h=(1/T)*mod(t,T);
u=compare(D,h);

RLesr=1e-3;
RCesr=1e-3;
Vd=0.7;
Rdson=1e-3;
RT=RLesr+Rdson;
ko=1+RCesr/R;

a11=-(1/L)*(  RT+(RCesr/ko)*(1-u)  );
a12=-(1/L)*(  ( 1-RCesr/(ko*R) )*(1-u)  );
a21= (1/C)*(  (1-u)/ko  );
a22=-(1/C)*(  1/(ko*R)  );
b1=  (Vin/L);

A=[a11 a12;a21 a22];
B=[b1;0];
dy=A*y+B;
end

function y = sat3(x)
y=1/2*(1+abs(x)-abs(x-1));
end

function a= compare(d,h)
b=(tanh((h-d)*1e12)+1)/2;
a=2-((tanh((b)*1e12)+1)/2)*2;
end