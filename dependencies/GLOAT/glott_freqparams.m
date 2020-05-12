function q=glott_freqparams(g,f0)
% GLOTTALFREQPARAMS Frequency-based parameterisation of the glottal flow
% Q=GLOTTALFREQPARAMS(G,F0) returns a structure containing multiple
% frequency-based parameters of the glottal flow.
%
% Description of the fields:
%
% DH12
%
% delta H_12 is more commonly known as H1-H2. It is the difference
% of the first two formants in decibels. See:
%
% I. Titze and J. Sundberg. Vocal intensity in speakers and
% singers. Journal of the Acoustical Society of America,
% 91(5):2936-2946, May 1992.
%
% PSP
%
% The parabolic spectrum parameter matches a second order
% polynomial to the flow spectrum computed over a single glottal
% cycle. See:
%
% P. Alku, H. Strik, and E. Vilkman. Parabolic spectral parameter
% -- a new method for quantifiction of the glottal flow. Speech
% Communication, 22:67-79.
%
% HRF
%
% Harmonic richness factor is the ratio of higher harmonics to the
% first harmonic. For more information, see:
%
% D. G. Childers and C. K. Lee. Vocal quality factors: Analysis,
% synthesis, and perception. Journal of the Acoustical Society of
% America, 90(5): 2394-2410.


% $Id: glottalfreqparams.m 127 2006-02-17 08:40:08Z mairas $

if nargin<2
  f0 = find_f0(g);
end

% get the fft of g

%G = fft(g);

DH12 = dh12(g,f0);

PSP = psp(g,f0);

HRF = hrf(g,f0);

q.DH12 = DH12;
q.PSP = PSP;
q.HRF = HRF;



%%%%

function y=dh12(g,f0)

Nfft = 2*1024;
h = hamming(length(g.s));
fft_dglot = 20*log10(abs(fft(h'.*g.s,Nfft)));

Nf0 = round(g.fs/f0);

F0 = round(Nfft/Nf0);
if isinf(F0)
    F0 = round(Nfft/100);
end


delta = fix(0.25*F0);
[h1 ind1] = max(fft_dglot(F0+1-delta:F0+1+delta));
[h2 ind2] = max(fft_dglot(2*F0+1-delta:2*F0+1+delta));
    
h12 = h1-h2;

% H1 = 20*log10(abs(G.at(f0)));
% H2 = 20*log10(abs(G.at(2*f0)));

y = h12;



%%%%
function PSP=psp(g,f0)

%tp = glott_timeparams(g,f0);
tp = get_closure_instants(g,f0);
T0 = 1/f0;

n0 = round(g.fs/f0);
%gmax = signal_ap(1/n0*ones(1,n0),g.fs);
gmax = 1/n0*ones(1,n0); %% added


gmax_zp = gmax; gmax_zp(2048)=0;
ws=warning('off','MATLAB:log:logOfZero');
%Xmax = 20*log10(half1(abs(fft(gmax_zp))));
Nfft = 2*1024; %% added
h = hamming(length(gmax_zp)); %% added
fft_gmax_zp = abs(fft(gmax_zp,Nfft)).^2; %% added

Xmax = 20*log10(fft_gmax_zp(1:Nfft/2)); % half1(abs(fft(gmax_zp))));
warning(ws);
amax = pspa_ap(Xmax);

PSP = [];
j=1;
for i=1:length(tp)
    
    if (tp(i)-n0 < 1)
        continue;
    end
    
  g_p = g.s(tp(i)-n0:tp(i));  %trim(g,tp(i).t_c-T0,tp(i).t_c);
  g_p0 = g_p-min(g_p);
  E_gp = sum(g_p.^2);
  g_pe = g_p0./(sqrt(E_gp));
  g_zp = g_pe; g_zp(2048)=0;
  %X = 20*log10(half1(abs(fft(g_zp))));
  fft_g_zp = abs(fft(g_zp,Nfft)).^2; %% added
  X = 20*log10(fft_g_zp(1:Nfft/2)); % half1(abs(fft(gmax_zp))));
  
  a = pspa_ap(X);
  
%   figure(1);
%   subplot(2,1,1); plot(g_zp);
%   subplot(2,1,2); plot(X);
%   waitforbuttonpress;
  
  PSP(j)=a/amax;
  j = j+1;
end

%%%%
function a=pspa_ap(X)

nel = 0.01;

done = 0;

N = 3;
while ~done
  k = 0:N-1;
  anum = N*sum(X(k+1).*k.^2)-sum(X(k+1)).*sum(k.^2);
  aden = N*sum(k.^4)-(sum(k.^2)).^2;
  a = anum / aden;
  
  b = 1/N*sum(X(k+1)-a*k.^2);
  
  NEnum = sum((X(k+1)-a*k.^2-b).^2);
  NEden = sum(X(k+1).^2);
  
  NE = NEnum/NEden;
  
  if NE<nel
    N = N+1;
  else
    done = 1;
  end
  
end

% 
%%%%

function tp = get_closure_instants(g,f0)

n0 = round(g.fs/f0);
T0 = 1/f0;

mx1 = max(g.s);
[pks, I_gmax] = findpeaks(g.s,'MinPeakDistance',n0*0.4,'MinPeakHeight',mx1*0.50);

T = I_gmax; %t_gmax;
tp = [];
j = 1;
for i=1:length(T)-1
    
    t_max=T(i);
    
    
    if (abs(length(g.s)-t_max)<1e-10)  || ((t_max+n0) > length(g.s))
        continue;
    end
    
    g_next = g.s(round(t_max):round((t_max+n0)));
    d_gnext=diff(g_next);
    [d_min,idmin]=min(d_gnext);
    t_dmin=(idmin);
    
    g_next2=g.s(t_dmin+t_max:t_max+n0);
    tm = find(diff(( diff(g_next2)  - 0)>0)>0)+1; 
    
    if length(tm)==0
        continue;
    end
    
    t_c = tm(1)+t_dmin+t_max-1;
    t_pb = t_c-n0;
    t_pe = t_c+round(0.2*n0);
    if (t_pb < 1) || (t_pe > length(g.s))
        continue;
    end
    
    tp(j) = t_c;
    j = j +1;
    
end


%%%%
function a=pspa(X)

nel = 0.01;

done = 0;

N = 3;
while ~done
  k = 0:N-1;
  anum = N*sum(X.s(k+1).*k.^2)-sum(X.s(k+1)).*sum(k.^2);
  aden = N*sum(k.^4)-(sum(k.^2)).^2;
  a = anum / aden;
  
  b = 1/N*sum(X.s(k+1)-a*k.^2);
  
  NEnum = sum((X.s(k+1)-a*k.^2-b).^2);
  NEden = sum(X.s(k+1).^2);
  
  NE = NEnum/NEden;
  
  if NE<nel
    N = N+1;
  else
    done = 1;
  end
end

%%%%
function HRF=hrf(g,f0)

 Nfft = 2*1024;
 h = hamming(length(g.s));
 %fft_dglot = (abs(fft(h.*dglot,Nfft))); % amplitude values
 fft_dglot = abs(fft(h'.*g.s,Nfft)).^2; % power spectrum % power spectrum is right one to be used
 
 Nf0 = round(g.fs/f0);
 
 F0 = round(Nfft/Nf0);
 %N = fix((Nfft/2)/(Nfft/Nf0));
 N = fix((Nfft*1500/g.fs)/(Nfft/Nf0));
    
 % Find peaks near the F0 values
 delta = fix(0.25*F0);
 h = zeros(N,1);
 ind = zeros(N,1);
 prevInd = 1;
 
 
 for i = 1:N
    [h(i) ind(i)] = max(fft_dglot(prevInd+F0-delta:prevInd+F0+delta));
    prevInd = prevInd+F0-delta+ind(i)-1; 
 %     plot(prevInd,fft_dglot(prevInd),'rx');
 end
    
    hr = 10*log10(sum(abs(h(2:end)))/abs(h(1))); 
    HRF = hr;

    
% Ga = abs(G);
% nh = get_harmonics(Ga,f0);
% 
% HRF = 20*log10(sum(Ga(nh(2:end)))/Ga.s(nh(1)));

% we can't calculate NRH unless we know the original window length

%%%%
function nh=get_harmonics(Xabs,f0)

% get N (vector of harmonics to get)
N = 1:floor((Xabs.frequency.fs/2)/f0);

nh=N;
for i=N
  % for each harmonic, find the local maximum in the region
  [f,n]=localmax(Xabs,f0,i);
  nh(i)=n;
end

%%%%
function [f,n]=localmax(Xabs,f0,N)
% LOCALMAX Find a largest spectral peak around N*f0
%   F=LOCALMAX(XABS,FX,F0,N) returns the frequency, at which XABS
%   gains its largest value. XABS is the absolute magnitude
%   spectrum, FX are the corresponding frequency points, F0 is the
%   assumed fundamental frequency and N is the order of harmonics.

fmin=(N-0.5)*f0;
fmax=(N+0.5)*f0;

Imin=at(Xabs,fmin);
Imax=at(Xabs,fmax);

[foo,I]=max(Xabs(Imin:Imax));

n=Imin-1+I;
f=Xabs.frequency.f(n);

