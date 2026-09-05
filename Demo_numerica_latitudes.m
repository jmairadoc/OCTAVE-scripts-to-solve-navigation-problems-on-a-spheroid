% Parámetros elipsoide
a = 6378137;
b = 6356752.314245;
k = b/a;

% Iteración a la milésima de grado (0.001°)

paso = 0.001;
lat_grados =  0 : paso : 90;
maxima_diferencia = 0;
latitud_maxima = 0;

% Bucle iterativo

for i = 1:length(lat_grados)
     phi_grados = lat_grados(i);

% Diferencia

    phi_radianes  = deg2rad(phi_grados);
    beta_radianes = atan(k * tan(phi_radianes));
    beta_grados = rad2deg(beta_radianes);
    dif_grados = phi_grados - beta_grados;

% Valor máximo

    if dif_grados >  maxima_diferencia

        maxima_diferencia = dif_grados;
        latitud_maxima= phi_grados;
    end

    disp ("grados: "), disp (phi_grados);
    disp ("diferencia: "), disp (dif_grados);

end
disp ("La diferencia máxima es:  "),  disp (maxima_diferencia);
disp ("en la latitud: "), disp  (latitud_maxima);
