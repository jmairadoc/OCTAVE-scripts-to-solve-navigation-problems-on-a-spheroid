%**************************************************************************************
%***********************          RESOLUCIÓN DEL          *****************************
%***********************  PROBLEMA DIRECTO DE NAVEGACION  *****************************
%***********************         EN UN ELIPSOIDE          *****************************
%**************************************************************************************

clear; clc;

%PREFERENCIAS OCTAVE
format long%formato con 15 dígitos significativos


%***********   PARÁMETROS DEL ELIPSOIDE (WGS84 a=6378137 b=6356752.314245)   **********
%a = input("Ingresa radio polar: ")
%b = input("Ingresa radio ecuatorial ")
a=6378137;%semieje mayor, radio ecuatorial
b=6356752.314245;%semieje menor, radio polar
f=(a-b)/a;%achatamiento
e=sqrt(f*(2-f));%primera excentricidad
eprima=e/sqrt(1-e^2);%segunda excentricidad
ecuadrado=f*(2-f);%primera excentricidad al cuadrado
eprimacuadrado=e^2/(1-e^2);%segunda excentricidad al cuadrado
c=sqrt((a^2/2)+(b^2/2)*(1/(tanh(e)*e)));%radio de curvatura polar
n=(a-b)/(a+b);%tercer achatamiento
%**************************************************************************************


%******************************** Datos de entrada ************************************
phi1       = input("Ingresa latitud de salida (S(-)N(+)): ");
lambda1    = input("Ingresa longitud de salida (W(-)E(+)): ");
alpha1     = input("Ingresa rumbo inicial (circular): ");
s12millas  = input("Ingresa distancia (millas): ");

%phi1=-10%latitud inicio dato S(-) N(+)
%lambda1=0%longitud inicio dato W(-) E(+)
%alpha1=210%azimuth inicio dato en circular
%s12millas=600%distancia dato en millas
%**************************************************************************************


%********************************** Millas a metros ***********************************
s12 = s12millas*1852;
%**************************************************************************************


%************************************ Cálculos ****************************************
beta1=atand((1-f)*tand(phi1));% beta1 es latitud reducida de phi1 en esfera auxiliar.
alpha0=asind(sind(alpha1)*cosd(beta1));%(5)

if phi1==0 %En este caso el triángulo NEA no existe, entonces:
  alpha0=alpha1;
  sigma1=0;
  w1=0;
else
  [sigma1, w1, alpha0] = funcion_NEA(beta1, alpha1); % funcion que resuelve NEA
endif

k=eprima*cosd(alpha0);%(9)
epsilon=(sqrt(1+k^2)-1)/(sqrt(1+k^2)+1);%(16)
A1=1/(1-epsilon)*(1+1/4*epsilon^2+1/64*epsilon^4+1/256*epsilon^6);%(17)

%hallando I1sigma1 ecuación(15)
[I1sigma1] = funcion_I1sigma1(epsilon, sigma1);

s1=b*I1sigma1;%(7)
s2=s1+s12;
t2=(s2/(b*A1))*180/pi;


%hallando sigma2 ecuación(20)
[sigma2] = funcion_sigma2(epsilon, t2);

% Triángulo esférico rectilátero NEB por Neper
% Datos NEB a=90 b=90-beta2 B=alpha0 c=sigma2
beta2=asind(cosd(alpha0)*sind(sigma2));%cos(90-beta2))=sen(90-alpha0)*sen(sigma2)

if beta2==0 %NEB será un triángulo esférico birrectángulo sigma2=w2
  w2=sign(sind(alpha1))*sigma2; %necesario ponerle signo al w2
else
  w2=atan2d(tand(alpha0)*tand(beta2),cosd(sigma2)/cosd(beta2));
end

alpha2=acosd(cosd(alpha0)*cosd(w2));

if (alpha1 > 180)
   alpha2=360-alpha2;
end

%hallando I3sigma1 e I3sigma2 ecuación(23)
[I3sigma1, I3sigma2] = funcion_I3sigma1_I3sigma2(epsilon, sigma1, sigma2, n);

%hallando lambda12 diferencia de longitud
lambda1g=(w1*pi/180-f*sind(alpha0)*(I3sigma1))*180/pi;%(8)
lambda2g=(w2*pi/180-f*sind(alpha0)*(I3sigma2))*180/pi;%(8)
lambda12=lambda2g-lambda1g;

lambda2=lambda1+lambda12;

if lambda2>180
          lambda2=lambda2-360;
elseif lambda2<-180
          lambda2=360+lambda2;
end

phi2=(atan(tand(beta2)/(1-f)))*180/pi;%(6)
%**************************************************************************************


%********************* Display cálculos de salida *************************************
disp ("latitud final: "), disp (phi2);
disp ("longitud final: "), disp (lambda2);
disp ("rumbo llegada (circular): "), disp (alpha2);
%**************************************************************************************
