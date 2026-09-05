%hallando I3sigma12

function [I3sigma12] = funcion_I3sigma12(epsilon, sigma12, A3, n)

        C31=(1/4-1/4*n)*epsilon+(1/8-1/8*n^2)*epsilon^2+(3/64+3/64*n-1/64*n^2)*epsilon^3+(5/128+1/64*n)*epsilon^4+3/128*epsilon^5;
        C32=(1/16-3/32*n+1/32*n^2)*epsilon^2+(3/64-1/32*n-3/64*n^2)*epsilon^3+(3/128+1/128*n)*epsilon^4+5/256*epsilon^5;
        C33=(5/192-3/64*n+5/192*n^2)*epsilon^3+(3/128-5/192*n)*epsilon^4+7/512*epsilon^5;
        C34=(7/512-7/256*n)*epsilon^4+7/512*epsilon^5;
        C35=21/2560*epsilon^5;
        C31bis=C31*sind(2*1*sigma12);
        C32bis=C32*sind(2*2*sigma12);
        C33bis=C33*sind(2*3*sigma12);
        C34bis=C34*sind(2*4*sigma12);
        C35bis=C35*sind(2*5*sigma12);
        sumac=C31bis+C32bis+C33bis+C34bis+C35bis;
        I3sigma12=A3*(sigma12*pi/180+sumac);%(23)

end
