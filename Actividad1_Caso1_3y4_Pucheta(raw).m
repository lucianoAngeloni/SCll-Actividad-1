clear; close all; clc;
mediciones=xlsread('Curvas_Medidas_RLC.xls');
figure(1);
plot(mediciones(:,1),mediciones(:,2)); %corriente
figure(2);
plot(mediciones(:,1),mediciones(:,3)); %vc

%Prueba
%  R=2.07;
%  L=0.000010;
%  C=0.00021757;
% 
% %Matrices
% A=[[-R/L -1/L]; [1/C 0]];
% B=[[1/L]; [0]];
% Cm=[1 0];
% D=[0];
% 
% sys1=ss(A,B,Cm,D);%sys1:modelo de espacio de estados con salida de corriente
% Cm=[0 1];
% sys2=ss(A,B,Cm,D);%sys2:modelo de espacio de estados con salida de voltaje en capacitor


%definicion de la entrada
u=zeros(1,1000);%devuelve una matriz de 1 por 1000
i=1;%iteracion
h=1e-4;%paso
t=0:h:(0.1-h);
vin=zeros(1,1000); %la entrada es siempre el escalon de 12v
 % p=0;
  signo=true;
for(i=100:1:1000)
    if mod(i,500)==0
        signo=not(signo);
    end
    if signo==1
        u(1,i)=12;
    end
    if signo==0
        u(1,i)=-12;
    end
     
end
%while(i<=h+1) 
 %   u(i)=vin; %la entrada es siempre el escalon de 12v
  %  p=p+0.04;  
   %if(p>1e-3)
    %   vin= vin*(-1);
     %  p=0;
   %end 
   %i=i+1;
%end 
figure(3);
plot(t,u),hold on;

%obtencion de RLC
%Metodo de chen
%kn = y(tn)/K-1 con K=12
k1=5.20835075/12-1; %t1=0.013  y(t1)=5,20835075
k2=8.602641821/12-1; %2t1=0.016 y(2t1)=8,602641821
k3=10.3015648/12-1; %3t1=0.019 y(3t1)=10,3015648
be=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;
alfa1=(k1*k2+k3-sqrt(be))/(2*(k1^2+k2));
alfa2=(k1*k2+k3+sqrt(be))/(2*(k1^2+k2));
beta=(k1+alfa2)/(alfa1-alfa2);
T1=-0.003/log(alfa1);
T2=-0.003/log(alfa2);
T3=beta*(T1-T2)+T1;

s=tf('s');
%G=12*(T3*s+1)/((T1*s+1)*(T2*s+1));
G=12*(T3*s+1)/((T1*s+1)*(T2*s+1));
G1=G/12;
[num,den]=tfdata(G1,'v');%devuelve el num y dem de la funcion de transf

figure(4);
[yaprox,taprox]=lsim(G,u/12,t);
plot(taprox,yaprox);

figure(5);
plot(taprox,yaprox), hold on;
plot(mediciones(:,1),mediciones(:,3)), hold on;

%DEDUCCION DE RLC
 [num,den] =     tfdata(G,'v');      %LA 'v' ES PARA QUE SE ME LO GUARDE EN UN VECTOR
 den_norm  =     den/den(1);         %DIVIDO TODO EL DENOMINADOR PARA NORMALIZAR EL COEFICIENTE DE S^2
 num_norm  =     num/(12*den(1));%NORMALIZO EL DENOMINADOR PARA ESCALON
 L         =     0.1;               %ASUMO UN VALOR DE L
 R         =     L*den_norm(2);
 C         =     1/(L*den_norm(3));

% R=2.07;
% L=0.000010;
% C=0.00021757;
% num2=s/L;
% %den2=((s^2)+(s*(R/(L)))+(1/(L*C)));
% den3=s*(R/(L));
% den4=(1/(L))*(1/(C));
% G_i=num2/den2;
G_i     =       (s/(L))/(s^2+(R/L)*s+(1/(L*C)));

figure(6);
plot(mediciones(:,1),mediciones(:,2)); hold on;
lsim(G_i,u,t, 'R');
axis([0 0.12 -0.1 0.06]);


%Verificacion de salida
%t1=[1:1:1000];
%t1=[];
%j=[500:1:1000];
%X0=[0 0]';x=[0 0]';
%h1=7.5151e-7;
%tf=0.05e-3;
%pasos=tf/h1;
%j=1;
%while(j<=pasos+1)
 %   t1(j)=j*h1;
  %  vin(j)=-12;
%     i1(j)=x(1);  
%     vc(j)=x(2);
%     xp=A*(x-X0)+B*vin(j); %xpunto 
%     x=x+(h1*xp); %para obtener x integro por Euler xp 
%     y=C*x; %salida
%     j=j+1;
% end
% 
% figure(6);
% plot(t1,i1);





