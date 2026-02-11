%hallando sigma2 ecuacion(20)
function [sigma2] = funcion_sigma2(epsilon, t2)

            Cprima11=1/2*epsilon-9/32*epsilon^3+205/1536*epsilon^5;
            Cprima12=5/16*epsilon^2-37/96*epsilon^4+1335/4096*epsilon^6;
            Cprima13=29/96*epsilon^3-75/128*epsilon^5;
            Cprima14=539/1536*epsilon^4-2391/2560*epsilon^6;
            Cprima15=3467/7680*epsilon^5;
            Cprima16=38081/61440*epsilon^6;
            Cpr11bis=Cprima11*sind(2*1*t2);
            Cpr12bis=Cprima12*sind(2*2*t2);
            Cpr13bis=Cprima13*sind(2*3*t2);
            Cpr14bis=Cprima14*sind(2*4*t2);
            Cpr15bis=Cprima15*sind(2*5*t2);
            Cpr16bis=Cprima16*sind(2*6*t2);
            sumab=Cpr11bis+Cpr12bis+Cpr13bis+Cpr14bis+Cpr15bis+Cpr16bis;
            sigma2=t2+sumab*180/pi;%(20)

end
