% Integral I1sigma2
function [I1sigma2] = funcion_I1sigma2(epsilon, sigma2)

          A1=1/(1-epsilon)*(1+1/4*epsilon^2+1/64*epsilon^4+1/256*epsilon^6);%17
          C11=-1/2*epsilon+3/16*epsilon^3-1/32*epsilon^5;
          C12=-1/16*epsilon^2+1/32*epsilon^4-9/2048*epsilon^6;
          C13=-1/48*epsilon^3+3/256*epsilon^5;
          C14=-5/512*epsilon^4+3/512*epsilon^6;
          C15=-7/1280*epsilon^5;
          C16=-7/2048*epsilon^6;
          C11bis=C11*sind(2*1*sigma2);
          C12bis=C12*sind(2*2*sigma2);
          C13bis=C13*sind(2*3*sigma2);
          C14bis=C14*sind(2*4*sigma2);
          C15bis=C15*sind(2*5*sigma2);
          C16bis=C16*sind(2*6*sigma2);
          suma=C11bis+C12bis+C13bis+C14bis+C15bis+C16bis;
          I1sigma2=A1*(sigma2*pi/180+suma);%(15)

end
