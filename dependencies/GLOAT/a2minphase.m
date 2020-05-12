function [ a ] = a2minphase( a, Nfft )
%A2MINPHASE Summary of this function goes here
%   Detailed explanation goes here
    
    p = length(a)-1;
    afft = abs(fft(a,Nfft));
    thresh = 0.000001;
   
    afft(afft < thresh) = thresh;

    afft = 1./afft;
    afft = afft.^2;
    ac = real(ifft(afft));
    a = levinson(ac,p);
    
end

