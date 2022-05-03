mediciones=xlsread('Curvas_Medidas_RLC.xls');
figure(1);
plot(mediciones(:,1),mediciones(:,2)); %corriente
figure(2);
plot(mediciones(:,1),mediciones(:,3)); %vc

L=1;
C=1;
R=1;

%Matrices
A=[[-R/L -1/L]; [1/C 0]];
B=[[1/L]; [0]];
C=[1 0];
D=[0];

sys1=ss(A,B,C,D);%sys1:modelo de espacio de estados con salida de corriente
C=[0 1];
sys2=ss(A,B,C,D);%sys2:modelo de espacio de estados con salida de voltaje en capacitor


%definicion de la entrada
i=1;%iteracion
h=1e-4;%paso
t=0:h:(0.1-h);
vin=12; %la entrada es siempre el escalon de 12v
  p=0;
while(i<=h+1) 
    u(i)=vin; %la entrada es siempre el escalon de 12v
    p=p+0.04;  
   if(p>1e-3)
       vin= vin*(-1);
       p=0;
   end 
   i=i+1;
end 
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


