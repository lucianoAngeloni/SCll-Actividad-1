clear; close all; clc;

%constantes punto 1
%R=4700;
%L=0.00001;
%C=0.0000001;
%punto2
R=5600;
L=0.00001;
C=0.0000001;

i=1;%iteracion
t=[];%arreglo

%se uso para calcula h y tf
A=[[-R/L -1/L]; [1/C 0]];
B=[[1/L]; [0]];
C=[R 0];
D=[0];
[num,den]=ss2tf(A,B,C,D,1);% eig(A) se puede hacer asi mas facil
roots(den);

%calculo de las constantes de tiempo
h=3.6e-11; % periodo de muestreo
tf=4.22e-3; % intervalo de integracion
pasos=tf/h; %
X0=[0 0]';x=[0 0]';  
  vin=12; %la entrada es siempre el escalon de 12v
  p = 0;
while(i<=pasos+1) 
    t(i)=i*h; %t aumenta con i 
    u(i)=vin; %la entrada es siempre el escalon de 12v
    p = p + h; 
   %variables de estado del sistema  
   i1(i)=x(1);  
   vc(i)=x(2);
   if(p>1e-3)
       vin= vin*(-1);
       p = 0;
   end
       
    
   %Sistema modelado en el espacio de estados 
   xp=A*(x-X0)+B*u(i); %xpunto 
    
   x=x+(h*xp); %para obtener x integro por Euler xp 
    
   y=C*x; %salida 
    
   i=i+1; 
end 
%Grafico 
subplot(3,1,1); plot(t,i1,'r'); title('I1,t');  
subplot(3,1,2); plot(t,vc,'b'); title('Vc,t'); 
subplot(3,1,3); plot(t,u,'g');  title('U1,t');






