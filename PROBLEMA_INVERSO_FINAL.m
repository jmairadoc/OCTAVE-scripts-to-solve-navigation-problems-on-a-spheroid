%************************************************************************************
%******************************          RESOLUCIÓN DEL          ********************
%******************************  PROBLEMA INVERSO DE NAVEGACION  ********************
%******************************         EN UN ELIPSOIDE          ********************
%************************************************************************************


%PREFERENCIAS OCTAVE
format long %formato con 15 dígitos significativos


%********   PARÁMETROS DEL ELIPSOIDE (WGS84 a=6378137 b=6356752.314245)   ***********
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
%************************************************************************************


%*************************************** Datos de entrada ***************************
%DATOS
phi1       = input("Ingresa latitud de salida: ");
lambda0    = input("Ingresa longitud de salida: ");
phi2       = input("Ingresa latitud final: ");
lambda2    = input("Ingresa longitud final: ");

%phi1=70;    %latitud inicial dato S(-) N(+)
%lambda0=10;   %longitud de salida
%phi2=-30;   %latitud final dato S(-) N(+)
%lambda2=-10;   %longitud de llegada (NO cerca del antipodal)
deltalambda=lambda2-lambda0;     %longitud de llegada menos longitud de salida
%************************************************************************************


%*********************************************** Simetrias **************************
simetria = 0;

if (phi1 == 0 && phi2 == 0) %Caso en el que se navega por el ecuador,(no es simetría)
  simetria = 1;

elseif (phi1 < 0 && phi2 < 0 && phi2 < phi1  && deltalambda < 0) %HS Caso rumbo SW
  [phi1, phi2] = deal(phi2, phi1);
  simetria = 2;

elseif (phi1 < 0 && phi2 < 0 && phi2 < phi1 && deltalambda >= 0) %HS Caso rumbo SE
  [phi1, phi2] = deal(phi2, phi1);
  simetria = 3;

elseif (phi1 >= 0 && phi2 >= phi1)   %HN Caso en el que se navega hacia NE o NW
  [phi1, phi2] = deal(-phi2, -phi1);
  simetria = 4;

elseif ((phi1 > 0 && phi2 > 0 && phi2 <= phi1 && deltalambda >= 0)  %HN Caso rumbo SE
                                  ||
       (phi1 > 0 && phi2 <= 0 && deltalambda >= 0))    %HN-HS Caso rumbo SE
  [phi1, phi2] = deal(-phi1, -phi2)
  simetria = 5;

elseif ((phi1 > 0 && phi2 > 0 && phi2 <= phi1  && deltalambda < 0)  %HN Caso rumbo SW
                                  ||
       (phi1 > 0 && phi2 <= 0  && deltalambda < 0))    %HN-HS Caso rumbo SW
  [phi1, phi2] = deal(-phi1, -phi2);
  simetria = 6;

endif
simetria;
%************************************************************************************


%**************** Primera aproximación para el método de Newton *********************
beta1=atand((1-f)*tand(phi1));%formula 6 Clairaut
beta2=atand((1-f)*tand(phi2));%formula 6 Clairaut

%Para el rumbo inicial se asume que es w12 = lambda12/wmedia
wmedia=sqrt(1-e^2*((cosd(beta1)+cosd(beta2))/2)^2);%48
w12=deltalambda/wmedia

%Resolución del triángulo esférico NAB
sigma12=acosd(sind(beta1)*sind(beta2)+cosd(beta1)*cosd(beta2)*cosd(w12))
alpha1=atan2d(sind(w12)*cosd(beta2),(cosd(beta1)*sind(beta2)-sind(beta1)*cosd(beta2)*cosd(w12)))
alpha2=atan2d(sind(w12)*cosd(beta1),(-cosd(beta2)*sind(beta1)+sind(beta2)*cosd(beta1)*cosd(w12)))
%Distancia recorrida en el esferoide
s12=a*wmedia*sigma12*(pi/180);
s12 = abs(s12);
%************************************************************************************


%*************************** Comienzo del bucle iterativo ***************************
for (i=0:1:20)

alpha0=asind(sind(alpha1)*cosd(beta1));%formula 5
k=eprima*cosd(alpha0);%9
epsilon=(sqrt(1+k^2)-1)/(sqrt(1+k^2)+1);%16

if (simetria == 1) % Break navegación por el ecuador.
    break
endif

% Triángulo esférico rectilátero NEA por Neper
% Datos NEA a=90 b=90-beta1 B=alpha0 A=180-alpha1
[sigma1, w1] = funcion_NEA(beta1, alpha0, alpha1, phi1);% Llama a la función_NEA
sigma1;
w1;
% Triángulo esférico rectilátero NEB por Neper
% Datos NEB a=90 b=90-beta2 B=alpha0  (sugerencia de usar la misma función NEA)
% EN EL INVERSO, TENEMOS COMO DATO PHI2 POR LO QUE SE APLICA ECUACIÓN (6) CLAIRAUT
% tanB2 = (1-f)tan(phi2)
w2=asind(tand(beta2)*tand(alpha0)); %utilizar el arcotan igualqueen NEA envezarcoseno
alpha2=atan2d(sind(alpha0)/cosd(beta2),cosd(alpha0)*cosd(w2)); %ecuación 45 karney
%cosalpha2=cosd(sqrt((cosd(alpha1))^2*(cosd(beta1))^2+((cos(beta2))^2-(cosd(beta1))^2))/cosd(beta2)) %(45) Karney
%alpha2=acosd(cosalpha2)
sigma2=atan2d(sind(beta2)/cosd(alpha0),cosd(beta2)*cosd(w2));

%hallando I3sigma1 e I3sigma2 ecuación(23)
[I3sigma1, I3sigma2] = funcion_I3sigma1_I3sigma2(epsilon, sigma1, sigma2, n);

%hallando lambda12 diferencia de longitud (calculada según alpha1)
lambda1g=(w1*pi/180-f*sind(alpha0)*(I3sigma1))*180/pi;%(8)
lambda2g=(w2*pi/180-f*sind(alpha0)*(I3sigma2))*180/pi;%(8)
lambda12=lambda2g-lambda1g;

%sigmalambda12
sigmalambda12 = lambda12 - deltalambda;
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
m12 = b*(sqrt(1+k^2*(sind(sigma2))^2)*cosd(sigma1)*sind(sigma2)-sqrt(1+k^2*(sind(sigma1))^2)*sind(sigma1)*cosd(sigma2)-cosd(sigma1)*cosd(sigma2)*(Jsigma2-Jsigma1));

%hallando sigmalambda12
sigmalambda12 = lambda12 - deltalambda;
error = abs(sigmalambda12);

%hallando ddelta12/dalpha1
dlambda12_dalpha1 = m12/(a*cosd(alpha2)*cosd(beta2));

if error < 0.000000009
  break
endif

%hallando deltagriegaalpha1
sigmaalpha1 = (-sigmalambda12)/(dlambda12_dalpha1);

%hallando nuevo alpha1
alpha1 = alpha1 + sigmaalpha1;
end
%********************************** fin bucle iterativo *****************************


%**************************** Actualizar s12 y condición ecuador ********************
if (simetria == 1)
  error = 0;
else
  s1 = b*I1sigma1;%(7)
  s2 = b*I1sigma2;%(7)
  s12 = s2-s1;
  s12 = abs(s12);
endif
%************************************************************************************


% *********************************** Deshacer simetría *****************************
if  (simetria == 2)
    [phi1, phi2] = deal(phi2, phi1);
    [alpha1, alpha2] = deal(-180-alpha2, -180-alpha1);

  elseif (simetria == 3)
    [phi1, phi2] = deal(phi2, phi1);
    [alpha1, alpha2] = deal(180-alpha2, 180-alpha1);

 elseif (simetria == 4)
    [phi1, phi2] = deal(-phi2, -phi1);
    [alpha1, alpha2] = deal(alpha2, alpha1);

  elseif (simetria == 5)
    [phi1, phi2] = deal(-phi1, -phi2);
    [alpha1, alpha2] = deal(180-alpha1, 180-alpha2);

  elseif (simetria == 6)
    [phi1, phi2] = deal(-phi1, -phi2);
    [alpha1, alpha2] = deal(-180-alpha1, -180-alpha2);
endif
simetria;
%************************************************************************************


%********************* Paso a rumbos circulares y a millas **************************
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
s12millas = (s12/1852); % Convierte metros en millas
%************************************************************************************


%**************************** Display cálculos de salida ****************************
disp ("latitud inicial: "), disp (phi1)
disp ("longitud inicial: "), disp (lambda0)
disp ("latitud final: "), disp (phi2)
disp ("longitud final: "), disp (lambda2)
disp ("rumbo salida: "), disp (alpha1)
disp ("en circular:"), disp(alpha1cir)
disp ("rumbo llegada: "), disp (alpha2)
disp ("en circular:"), disp(alpha2cir)
disp ("distancia recorrida: "), disp (s12)
disp ("en millas: "), disp (s12millas)
disp ("error: "), disp (error)
disp ("iteraciones: "), disp(i)
%************************************************************************************
