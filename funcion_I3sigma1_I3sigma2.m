%hallando I3sigma1 e I3sigma2 ecuación(23)

function [I3sigma1, I3sigma2] = funcion_I3sigma1_I3sigma2(epsilon, sigma1, sigma2, n)
        A3=1-(1/2-1/2*n)*epsilon-(1/4+1/8*n-3/8*n^2)*epsilon^2-(1/16+3/16*n+1/16*n^2)*epsilon^3-(3/64+1/32*n)*epsilon^4-3/128*epsilon^5;%(24);
        C31=(1/4-1/4*n)*epsilon+(1/8-1/8*n^2)*epsilon^2+(3/64+3/64*n-1/64*n^2)*epsilon^3+(5/128+1/64*n)*epsilon^4+3/128*epsilon^5;
        C32=(1/16-3/32*n+1/32*n^2)*epsilon^2+(3/64-1/32*n-3/64*n^2)*epsilon^3+(3/128+1/128*n)*epsilon^4+5/256*epsilon^5;
        C33=(5/192-3/64*n+5/192*n^2)*epsilon^3+(3/128-5/192*n)*epsilon^4+7/512*epsilon^5;
        C34=(7/512-7/256*n)*epsilon^4+7/512*epsilon^5;
        C35=21/2560*epsilon^5;
        C31bis=C31*sind(2*1*sigma1);
        C32bis=C32*sind(2*2*sigma1);
        C33bis=C33*sind(2*3*sigma1);
        C34bis=C34*sind(2*4*sigma1);
        C35bis=C35*sind(2*5*sigma1);
        sumac=C31bis+C32bis+C33bis+C34bis+C35bis;
        I3sigma1=A3*(sigma1*pi/180+sumac);%(23)
        C31dis=C31*sind(2*1*sigma2);
        C32dis=C32*sind(2*2*sigma2);
        C33dis=C33*sind(2*3*sigma2);
        C34dis=C34*sind(2*4*sigma2);
        C35dis=C35*sind(2*5*sigma2);
        sumad=C31dis+C32dis+C33dis+C34dis+C35dis;
        I3sigma2=A3*(sigma2*pi/180+sumad);%(23)

end
