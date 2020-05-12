function [ gci_ins, goi_ins ] = gci( x, meanf0, fs )
%  Estimates the GCI (and GOI) instants of the given frame x according to
%  the GLOAT algorithm (Drugman, T and Dutoit, T: "Glottal Closure and
%  Opening Instant Detection from Speech Signals", 2009)
%  meanf0 = average f0 value of the speech sample (in Hz)
%  fs = sampling frequency (in Hz)



N = length(x);


 % Window length for mean-based signal
if 1/meanf0 < 0.008
    Tf0 = fs/meanf0;
 Nw =round(1.4*Tf0/2);
else
    Nw = 0.006*fs;
end
x = x(:)';

% LPC order
p = 10; 
if fs == 16000
  p = 18;  
elseif fs == 48000;
    p = 72;    
end
%p = round(fs/1000+2);



% Compute mean-based signal

x2 = [zeros(1,Nw), x, zeros(1,Nw+1)];
w = blackman(2*Nw+1);

y = zeros(1,N);

for n = 1:N
    for m=-Nw:Nw
   y(n) = y(n) + w(m+Nw+1)*x2(n+m+Nw+1);    
    end
    y(n) = y(n)/(2*Nw+1);
end

% Compute the derivative of the mean-based signal
ydiff = gradient(y);

% Compute GCI and GOI intervals, format: columns = one interval, rows =
% interval start and end samples

goi_int = [];
gci_int = [];

for i = 1:N-1
   if sign(ydiff(i)*ydiff(i+1)) == -1
       if sign(ydiff(i)) == 1
           for j=i:N-1
              if sign(y(j)*y(j+1)) == -1
                  break
              end
           end
           goi_int = [goi_int [i; j]];
       else
           for j=i:N-1
              if sign(y(j)*y(j+1)) == -1
                  break
              end
           end
           gci_int = [gci_int [i; j]];
       end
   end
end

% Find maximum values of the LPC residual from the intervals

x_res = filter(lpc(x,p),1,x);

gci_ins = zeros(1,size(gci_int,2));
goi_ins = zeros(1,size(goi_int,2));

for i=1:size(gci_int,2)
   [peak ind] = max(x_res(gci_int(1,i):gci_int(2,i)));
   gci_ins(i) = gci_int(1,i)+ind-1;
end


for i=1:size(goi_int,2)
   [peak ind] = max(x_res(goi_int(1,i):goi_int(2,i)));
   goi_ins(i) = goi_int(1,i)+ind-1;
end


% Plot results (uncomment below to plot)

% figure(2);
% 
% plot(x_res);
% 
% hold on;
% plot(x,'g');
% plot(100*y,'r');
% plot(goi_ins,x_res(goi_ins),'ro');
% plot(gci_ins,x_res(gci_ins),'rx');
% hold off;

end

