%***************************************************************************
%*********************          RESOLUCIÓN DEL          ********************
%*********************  PROBLEMA INVERSO DE NAVEGACION  ********************
%*********************         EN UN ELIPSOIDE          ********************
%***************************************************************************

clear; clc;

%PREFERENCIAS OCTAVE
format long %formato con 15 dígitos significativos


%*****   PARÁMETROS DEL ELIPSOIDE (WGS84 a=6378137 b=6356752.314245)   *****
%a = input("Ingresa radio polar: ")
%b = input("Ingresa radio ecuatorial: ")
a=6378137;%semieje mayor, radio ecuatorial
b=6356752.314245;%semieje menor, radio polar
f=(a-b)/a;%achatamiento
e=sqrt(f*(2-f));%primera excentricidad
eprima=e/sqrt(1-e^2);%segunda excentricidad
ecuadrado=f*(2-f);%primera excentricidad al cuadrado
eprimacuadrado=e^2/(1-e^2);%segunda excentricidad al cuadrado
c=sqrt((a^2/2)+(b^2/2)*(1/(tanh(e)*e)));%radio de curvatura polar
n=(a-b)/(a+b);%tercer achatamiento
%***************************************************************************



%***************************** Datos de entrada ****************************
%DATOS
phi1       = input("Ingresa latitud de salida: ");
lambda1    = input("Ingresa longitud de salida: ");
phi2       = input("Ingresa latitud final: ");
lambda2    = input("Ingresa longitud final: ");

%phi1=10;    %latitud inicial dato S(-) N(+)
%lambda1=0;   %longitud de salida
%phi2=-5;   %latitud final dato S(-) N(+)
%lambda2=-90;   %longitud de llegada (NO cerca del antipodal)
%***************************************************************************

%*************** Condición: diferencia de longitud entre 0 y pi ************
lambda12real = lambda2 - lambda1;
if (lambda12real > 180)
    lambda12real = -360 + lambda12real
elseif (lambda12real < -180)
        lambda12real = 360 - abs(lambda12real)
endif
%****************************************************************************

%******************************** Simetrias *********************************
simetria = 0; %No se aplican simetrías

if (phi1 == 0 && phi2 == 0) %Caso en el que se navega por el ecuador
  simetria = 1;

elseif ((phi1 <= 0 && phi2 < phi1) %simetría 2 y 7 (meridiano)
                   ||
       (phi1 > 0   && phi2 < 0 && phi2 < -phi1))
    [phi1, phi2] = deal(phi2, phi1);
    simetria = 2;

elseif ((phi1 <=0 && phi2 > -phi1) %simetría 3 y 5(meridiano + ecuador)
                 ||
       (phi1 > 0 && phi2 > 0 && phi2 > phi1))
    [phi1, phi2] = deal(-phi2, -phi1);
    simetria = 3;

elseif ((phi1 > 0 && phi2 > 0 && phi2 < phi1) %simetría 4 y 6 (ecuador)
                 ||
       (phi1 > 0 && phi2 < 0 && phi2 > -phi1))
    [phi1, phi2] = deal(-phi1, -phi2);
    simetria = 4;

endif

%***************************************************************************

%************ Primera aproximación para el método de Newton ****************
beta1=atand((1-f)*tand(phi1));%formula 6 Clairaut
beta2=atand((1-f)*tand(phi2));%formula 6 Clairaut

%Para el rumbo inicial se asume que es w12 = lambda12real/wmedia
wmedia=sqrt(1-e^2*((cosd(beta1)+cosd(beta2))/2)^2);%48
w12=lambda12real/wmedia;

%Resolución del triángulo esférico NAB
sigma12=acosd(sind(beta1)*sind(beta2)+cosd(beta1)*cosd(beta2)*cosd(w12));
alpha1=atan2d(sind(w12)*cosd(beta2),(cosd(beta1)*sind(beta2)
                                        -sind(beta1)*cosd(beta2)*cosd(w12)));
alpha2=atan2d(sind(w12)*cosd(beta1),(-cosd(beta2)*sind(beta1)
                                        +sind(beta2)*cosd(beta1)*cosd(w12)));
s12=a*wmedia*sigma12*(pi/180);
s12 = abs(s12);
%***************************************************************************

%******************** Comienzo del bucle iterativo *************************
for (i=0:1:20)

alpha0=asind(sind(alpha1)*cosd(beta1));%formula 5
k=eprima*cosd(alpha0);%9
epsilon=(sqrt(1+k^2)-1)/(sqrt(1+k^2)+1);%16

if (simetria == 1) % Break navegación por el ecuador.
    break
endif

% Triángulo esférico rectilátero NEA por Neper
% Datos NEA a=90 b=90-beta1 B=alpha0 A=180-alpha1
if alpha1==90 && phi1==0 %En este caso el triángulo NEA "no existe"
    alpha0=90;
    sigma1=0;
    w1=0;
  elseif alpha1==270 && phi1==0 %En este caso el triángulo NEA "no existe"
    alpha0=-90;
    sigma1=0;
    w1=0;
  else
    sigma1=atan2d((sind(beta1)/cosd(alpha0)),(cosd(beta1)*cosd(alpha1)/cosd(alpha0)));
    w1=atan2d(tand(beta1)*tand(alpha0),cosd(alpha1)/cosd(alpha0));
  endif
sigma1;
w1;
% Triángulo esférico rectilátero NEB por Neper
% Datos NEB a=90 b=90-beta2 B=alpha0
% EN EL INVERSO, DATO PHI2 POR LO QUE SE APLICA ECUACIÓN (6) CLAIRAUT
% tanB2 = (1-f)tan(phi2)
w2=asind(tand(beta2)*tand(alpha0));
alpha2=atan2d(sind(alpha0)/cosd(beta2),cosd(alpha0)*cosd(w2)); %en vez de(45)
sigma2=atan2d(sind(beta2)/cosd(alpha0),cosd(beta2)*cosd(w2));

%hallando I3sigma1 e I3sigma2 ecuación(23)
[I3sigma1, I3sigma2] = funcion_I3sigma1_I3sigma2(epsilon, sigma1, sigma2, n);

%hallando lambda12 diferencia de longitud (calculada según alpha1)
lambda1g=(w1*pi/180-f*sind(alpha0)*(I3sigma1))*180/pi;%(8)
lambda2g=(w2*pi/180-f*sind(alpha0)*(I3sigma2))*180/pi;%(8)
lambda12=lambda2g-lambda1g;

%I1Sigma1
[I1sigma1] = funcion_I1sigma1(epsilon, sigma1);
%I2Sigma1
[I2sigma1] = funcion_I2sigma1(epsilon, sigma1);
%Jsigma1
Jsigma1=I1sigma1-I2sigma1;
%I1Sigma2
[I1sigma2] = funcion_I1sigma2(epsilon, sigma2);
%I2Sigma2
[I2sigma2] = funcion_I2sigma2(epsilon, sigma2);
%Jsigma2
Jsigma2=I1sigma2-I2sigma2;
%m12
m12 = b*(sqrt(1+k^2*(sind(sigma2))^2)*cosd(sigma1)*sind(sigma2)
        -sqrt(1+k^2*(sind(sigma1))^2)*sind(sigma1)*cosd(sigma2)
        -cosd(sigma1)*cosd(sigma2)*(Jsigma2-Jsigma1));

%hallando deltalambda12
deltalambda12 = lambda12 - lambda12real;
error = abs(deltalambda12);

%hallando dlambda12/dalpha1
dlambda12_dalpha1 = m12/(a*cosd(alpha2)*cosd(beta2));

if error < 0.000000009
  break
endif

%hallando deltaaalpha1
deltaalpha1 = (-deltalambda12)/(dlambda12_dalpha1);

%hallando nuevo alpha1
alpha1 = alpha1 + deltaalpha1;
end
%*********************** fin bucle iterativo *******************************


%*************** Actualizar s12 y condición ecuador ************************
if (simetria == 1)
  error = 0;
else
  s1 = b*I1sigma1;%(7)
  s2 = b*I1sigma2;%(7)
  s12 = s2-s1;
  s12 = abs(s12);
endif
%***************************************************************************


%************************** Deshacer simetría ******************************
if  (simetria == 2)
    [phi2, phi1] = deal(phi1, phi2);
       [alpha1, alpha2] = deal(-180-alpha2, -180-alpha1);

 elseif (simetria == 3)
    [phi2, phi1] = deal(-phi1, -phi2);
    [alpha1, alpha2] = deal(alpha2, alpha1);

 elseif (simetria == 4)
    [phi1, phi2] = deal(-phi1, -phi2);
    [alpha1, alpha2] = deal(180-alpha1, 180-alpha2);

endif

%***************************************************************************


%***************** Paso a rumbos circulares y a millas *********************
if (alpha1 < 0)
     alpha1cir = (360 + alpha1);
  else
     alpha1cir = (alpha1);
endif

if (alpha2 < 0)
     alpha2cir = (360 + alpha2);
  else
     alpha2cir = (alpha2);
endif
s12millas = (s12/1852); % Convierte metros a millas
%***************************************************************************


%********************* Display cálculos de salida **************************
disp ("latitud inicial: "), disp (phi1)
disp ("longitud inicial: "), disp (lambda1)
disp ("latitud final: "), disp (phi2)
disp ("longitud final: "), disp (lambda2)
%disp ("rumbo salida: "), disp (alpha1)
disp ("rumbo salida (circular):"), disp(alpha1cir)
%disp ("rumbo llegada: "), disp (alpha2)
disp ("rumbo llegada (circular):"), disp(alpha2cir)
disp ("distancia recorrida (metros): "), disp (s12)
disp ("distancia recorrida (millas): "), disp (s12millas)
disp ("error: "), disp (error)
disp ("iteraciones: "), disp(i)
%***************************************************************************
