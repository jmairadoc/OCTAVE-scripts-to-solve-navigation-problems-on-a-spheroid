% ******************************************************
% **** Triángulo esférico rectilátero NEA por Neper ****
% **** Datos entrada: beta1, alpha1 ********************
% **** Datos salida: sigma1, w1, alpha0 ****************
% ******************************************************

function [sigma1, w1, alpha0] = funcion_NEA(beta1, alpha1);

alpha0=atan2d(sind(alpha1)*cosd(beta1),sign(beta1)*sqrt(1-(cosd(beta1)*sind(alpha1))^2));
sigma1=atan2d((sind(beta1)/cosd(alpha0)),(cosd(beta1)*cosd(alpha1)/cosd(alpha0)));
w1=atan2d(tand(beta1)*tand(alpha0),cosd(alpha1)/cosd(alpha0));

end
