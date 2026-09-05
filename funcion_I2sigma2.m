% Integral I2sigma2
function [I2sigma2] = funcion_I2sigma2(epsilon, sigma2)

            A2=(1-epsilon)*(1+1/4*epsilon^2+9/64*epsilon^4+25/256*epsilon^6);%42
            C21=1/2*epsilon+1/16*epsilon^3+1/32*epsilon^5;
            C22=3/16*epsilon^2+1/32*epsilon^4+35/2048*epsilon^6;
            C23=5/48*epsilon^3+5/256*epsilon^5;
            C24=35/512*epsilon^4+7/512*epsilon^6;
            C25=63/1280*epsilon^5;
            C26=77/2048*epsilon^6;
            C21bis=C21*sind(2*1*sigma2);
            C22bis=C22*sind(2*2*sigma2);
            C23bis=C23*sind(2*3*sigma2);
            C24bis=C24*sind(2*4*sigma2);
            C25bis=C25*sind(2*5*sigma2);
            C26bis=C26*sind(2*6*sigma2);
            suma=C21bis+C22bis+C23bis+C24bis+C25bis+C26bis;
            I2sigma2=A2*(sigma2*(pi/180)+suma);%(41)

end
