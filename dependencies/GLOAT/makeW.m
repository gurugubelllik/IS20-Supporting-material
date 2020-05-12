function [ w ] = makeW( x, p, DQ, PQ, d, Nramp, gci_ins, fs )
%   Create a AME weight function for frame x for LPC order p
%   DQ = duration quotient (from 0 to 1)
%   PQ = Position Quotient (from 0 to 1)
%   d = minimum value of the weight function
%   Nramp = length of the linear ramp (in samples)
%   gci_ins = Glottal Closure Instants of the frame x


N = length(x);
%Nramp = 1;
if Nramp > 0
UPramp = linspace(d,1,2+Nramp);
UPramp = UPramp(2:end-1);
DOWNramp = UPramp(end:-1:1);
end


if DQ+PQ > 1
    DQ = 1-PQ;
end

w = d.*ones(1,N+p);


try
for i = 1:length(gci_ins)-1
   T = gci_ins(i+1)-gci_ins(i);
   T1 = round(DQ*T);
   T2 = round(PQ*T);
   while T1+T2 > T
       T1 = T1-1;
   end
   w(gci_ins(i)+T2:gci_ins(i)+T2+T1-1) = 1;
   if Nramp > 0
       w(gci_ins(i)+T2:gci_ins(i)+T2+Nramp-1) = UPramp;
       if gci_ins(i)+T2+T1-Nramp > 0
           w(gci_ins(i)+T2+T1-Nramp:gci_ins(i)+T2+T1-1) = DOWNramp;
       end
   end
end

Nend = N-(T2+gci_ins(i+1));

if T2+gci_ins(i+1) < N
    if T1+T2 < Nend
        w(gci_ins(i+1)+T2:gci_ins(i+1)+T2+T1-1) = 1;
        if Nramp > 0
            w(gci_ins(i+1)+T2:gci_ins(i+1)+T2+Nramp-1) = UPramp;
            w(gci_ins(i+1)+T2+T1-Nramp:gci_ins(i+1)+T2+T1-1) = DOWNramp;
        end
    else
        T1 = Nend-T2;
                w(gci_ins(i+1)+T2:gci_ins(i+1)+T2+T1-1) = 1;
        if Nramp > 0
            w(gci_ins(i+1)+T2:gci_ins(i+1)+T2+Nramp-1) = UPramp;
        end
    end
end

    
end

end

