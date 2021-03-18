 function [Ipv,I0,Rs,Rsh] = param_1D_2R_L(Isc,Voc,Imp,Vmp,a)
 global Vt

    A=-(2*Vmp-Voc)/(a*Vt)+(Vmp*Isc-Voc*Imp)/(Vmp*Isc+Voc*(Imp-Isc))
    B=-Vmp*(2*Imp-Isc)/(Vmp*Isc+Voc*(Imp-Isc))
    C=a*Vt/Imp
    D=(Vmp-Voc)/(a*Vt)
    
    M1=0.3361;
    M2=-0.0042;
    M3=-0.0201;
    sigma = -1-log(-B)-A;
    Wn =-1-sigma -2/M1* (1-1./(1+M1*sqrt(sigma/2)./(1+M2*sigma.*exp(M3*sqrt(sigma)))) );
    
    Rs=C*(Wn-(D+A));
    Rsh=(Vmp-Imp*Rs)*(Vmp-Rs*(Isc-Imp)-a*Vt)/((Vmp-Imp*Rs)*(Isc-Imp)-a*Vt*Imp);
    Ipv=(Rsh+Rs)/Rsh*Isc;
    I0=((Rsh+Rs)/Rsh*Isc-Voc/Rsh)/(exp((Voc)/(a*Vt)));

    


