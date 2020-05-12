function I_feat=Intonation_features(data,fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code to find voice quality features from the speech Signal
% Authors:  G Krishna  (Speech Lab, IIIT Hyderabad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------       inputs         --------
% data      : input raw speech signal
% fs        : sampling freqnecy
% --------      outputs         --------
% I_feat    : voice quality features (76 dimentional)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dependencies
%       bufferaroundgci.m
%       zff_analysis.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% (1) http://www.public.asu.edu/~visar/IS2018Supp.pdf
%
% Useful references:
% 
% [1] A. Tsanas: "Accurate telemonitoring of Parkinson's disease symptom
%     severity using nonlinear speech signal processing and statistical
%     machine learning", D.Phil. thesis, University of Oxford, 2012
%
% [2] A. Tsanas, M.A. Little, P.E. McSharry, L.O. Ramig: "Nonlinear speech 
%     analysis algorithms mapped to a standard metric achieve clinically 
%     useful quantification of average Parkinson�s disease symptom severity", 
%     Journal of the Royal Society Interface, Vol. 8, pp. 842-855, June 2011
%
% [3] A. Tsanas, M.A. Little, P.E. McSharry, L.O. Ramig: "New nonlinear 
%     markers and insights into speech signal degradation for effective 
%     tracking of Parkinson�s disease symptom severity", International
%     Symposium on Nonlinear Theory and its Applications (NOLTA), 
%     pp. 457-460, Krakow, Poland, 5-8 September 2010
%
% [4] A. Tsanas, M.A. Little, P.E. McSharry, L.O. Ramig: �Enhanced classical
%     dysphonia measures and sparse regression for telemonitoring of 
%     Parkinson's disease progression�, IEEE Signal Processing Society, 
%     International Conference on Acoustics, Speech and Signal Processing 
%     (ICASSP), pp. 594-597, Dallas, Texas, US, 14-19 March 2010
%
% [5] A. Tsanas, M.A. Little, C. Fox, L.O. Ramig: "Objective automatic 
%     assessment of rehabilitative speech treatment in Parkinson�s disease",
%     IEEE Transactions on Neural Systems and Rehabilitation Engineering, 
%     Vol. 22, pp. 181-190, January 2014
%
% [6] A. Tsanas: "Automatic objective biomarkers of neurodegenerative 
%     disorders using nonlinear speech signal processing tools", 8th 
%     International Workshop on Models and Analysis of Vocal Emissions for 
%     Biomedical Applications (MAVEBA), pp. 37-40, Florence, Italy, 16-18 
%     December 2013 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pre-process data
data = data - mean(data);

s = data/max(abs(data));

[vnvsig,~,zf] = zff_analysis(s,fs);

[~,~,zf1,gci,es,f0] = zff_analysis(zf,fs);

vgci=gci(find(vnvsig(gci)==1)); vf0=f0(find(vnvsig(gci)==1));

Pitchperiod=round(8000/mean(vf0));

A0=max(bufferaroundgci(abs(s),gci,round(Pitchperiod/3),round(Pitchperiod/3))');

F0_jitter_feat = jitter_shimmer(f0);                % F0 Jitter measures

Soe_shimmer_feat = jitter_shimmer(es);              % SOE Shimmer measures

A0_shimmer_feat = jitter_shimmer(A0);               % SOE Shimmer measures

F0_feat=F0_measures(f0);                            % F0 stat. measures

[HNR,NHR] = HNRFun(s,fs);                           % Harmonics to Noise Ratio (HNR) and Noise to Harmonics Ratio (NHR)

 [VF_close, VF_open] = dypsa(s,fs);
% 
 [GQ] = glottis_quotient(VF_close,VF_open, fs, 50, 500, 1);   % Glottis open/closed quotient [Glottis Quotient (GQ) in the JRSI paper]
% 
 [GNE] = GNE_measure(s,fs);                       % Glottal to Noise Excitation (GNE) related measures

PPE=F0_dispersion(f0,fs)

I_feat=[F0_feat F0_jitter_feat Soe_shimmer_feat A0_shimmer_feat HNR NHR PPE GQ GNE];

end


function PPE=F0_dispersion(F0,fs)

% PPE2 (improved version compared to Max's)
F0mean_healthy_control = 120; % mean(F0); % strictly speaking, adjust this for healthy controls males = 120 and females = 190
logF0signal = log_bb(F0/F0mean_healthy_control); % log transform F0 is beneficial - see my ICASSP2010 paper, also controlling for F0 in healthy controls
ARcoef = arcov(logF0signal, 10); % identify AR coefficients
sig_filtered = filter(ARcoef, 1, logF0signal);
sig_filtered = sig_filtered(round(0.001*fs):end);
% Obtain measure of dispersion - use simple histogramic means for now
PPEd = hist(sig_filtered, linspace(min(sig_filtered), max(sig_filtered), 100));
PPE = H_entropy(PPEd/sum(PPEd))/log_bb(length(PPEd));


end

function F0_feat=F0_measures(F0)

F0_feat=[mean(F0) median(F0) std(F0) min(F0) max(F0)];

end

function [measures] = jitter_shimmer(A0)

mean_Ampl = mean(A0);

% Mean absolute difference of successive cycles
measures(1) = mean(abs(diff(A0)));

% Mean absolute difference of successive cycles - expressed in percent (%)
measures(2) = 100*mean(abs(diff(A0)))/mean_Ampl;

% Perturbation quotient
[Ampl_PQ3] = perq1(A0,3);
measures(3) = Ampl_PQ3.classical_Schoentgen;
measures(4) = Ampl_PQ3.classical_Baken;
measures(5) = Ampl_PQ3.generalised_Schoentgen;

[Ampl_PQ5]=perq1(A0,5); % Use 5 cycle samples (Schoentgen)
measures(6)=Ampl_PQ5.classical_Schoentgen;
measures(7)=Ampl_PQ5.classical_Baken;
measures(8)=Ampl_PQ5.generalised_Schoentgen;

[Ampl_PQ11]=perq1(A0,11); % Use 11 cycle samples (Schoentgen)
measures(9)=Ampl_PQ11.classical_Schoentgen;
measures(10)=Ampl_PQ11.classical_Baken;
measures(11)=Ampl_PQ11.generalised_Schoentgen;

% zeroth order perturbation
measures(12) = mean(abs(A0-mean_Ampl));

% Shimmer(dB)
measures(13) = mean(20*(abs(log10((A0(1:end-1))./(A0(2:end))))));

% CV
measures(14)=mean((diff(A0)).^2)/(mean_Ampl)^2;

% TKEO
measures(15) = mean(abs(TKEO(A0)));
measures(16) = std(TKEO(A0));
Ampl_TKEO_prc = prctile(TKEO(A0),[5 25 50 75 95]);
measures(17)=Ampl_TKEO_prc(1);
measures(18)=Ampl_TKEO_prc(2);
measures(19)=Ampl_TKEO_prc(3);
measures(20)=Ampl_TKEO_prc(4);

% AM
measures(21) = (max(A0)-min(A0))/(max(A0)+min(A0));

measures(22) = Ampl_TKEO_prc(4) - Ampl_TKEO_prc(1);
end


function PQ = perq1(time_series, K)

%% Calculate the PQ using the classical PQ formula

N = length(time_series);
mean_tseries = mean(time_series);
K1=round(K/2);
K2=K-K1;
p = 5;
sum1=0;

for i = K1:N-K2
    sum1 = sum1+mean(abs([time_series(i-K2:i+K2)]-time_series(i)));
end
        
PQ.classical_Schoentgen = (sum1/(N-K+1))/(mean_tseries);

sum2=0;
for i = K1:N-K2
    sum2 = sum2+mean(abs([time_series(i-K2:i+K2)]))-time_series(i);
end
        
PQ.classical_Baken = (sum2/(N-K+1))/(mean_tseries);

% perturbation quotient of the residue
time_series=time_series(:);
sum3=0;
% calculate the AR coefficients (I use the Yule-Walker equations)
new_tseries=(time_series-mean_tseries)';
a = aryule(time_series-mean_tseries,p);

for i = 1+p:N
    sum3 = sum3+abs(sum(a.*(new_tseries(i:-1:i-p))));
end

PQ.generalised_Schoentgen = (sum3/(N-p))/(mean_tseries);
    
end

function [energy] = TKEO(x)

data_length=length(x);
energy=zeros(data_length,1);

energy(1)=(x(1))^2; % first sample in the vector sequence

for n=2:data_length-1
    energy(n)=(x(n))^2-x(n-1)*x(n+1); % classical TKEO equation
end

energy(data_length)=(x(data_length))^2; % last sample in the vector sequence

end % end of TKEO function


function [GQ] = glottis_quotient(VF_close, VF_open, fs, f0min, f0max, flag)
%% Calculate the glottis quotients

cycle_open=abs(VF_open(2:end)-VF_close(1:end-1));
cycle_closed=abs(VF_open(1:end-1)-VF_close(1:end-1));

% remove erroneous cycles
if flag
    low_lim=fs/f0max; % lower limit
    up_lim=fs/f0min;  % upper limit
    N=length(cycle_open);
    for i=1:N-1
        if((cycle_open(i) > up_lim) || (cycle_open(i) < low_lim))
            cycle_open(i)=NaN;
        end
        if((cycle_closed(i) > up_lim) || (cycle_closed(i) < low_lim))
            cycle_closed(i)=NaN;
        end
    end
end

%statistics in time
prc1=prctile(cycle_open,[5 95]);
cycle_open_range_5_95_perc=prc1(2)-prc1(1);
prc2=prctile(cycle_closed,[5 95]);
cycle_closed_range_5_95_perc=prc2(2)-prc2(1);

GQ(1) = (cycle_open_range_5_95_perc/(cycle_open_range_5_95_perc+cycle_closed_range_5_95_perc));
GQ(2) = (nanstd(cycle_open));
GQ(3) = (nanstd(cycle_closed));

end


function [GNE] = GNE_measure(data, fs)

filt_order=100;
new_fs=10000;
x = 0.03*new_fs;
tstep=0.01*new_fs;
BW = 1000; % bandwidth
Fshift = 500; % shift fr   
data = resample(data, new_fs, fs); 

% cut-off1
Fc1=1:Fshift:(new_fs/2-BW-500); % not to cross Nq freq!
% cut-off2
Fc2=Fc1+BW;

for j=1:length(Fc1);
    d(j) =fdesign.bandpass('n,fc1,fc2',filt_order,Fc1(j),Fc2(j),new_fs);
    hd(j) = design(d(j));
end

steps=(length(data)-x)/tstep;

for i=1:steps+1   
    tseries = data(1+(i-1)*tstep:(i-1)*tstep+x);
    Dwindow = hann(length(tseries));
    segment_sig = tseries.*Dwindow;

    a = lpc(segment_sig,13);
    est_x = filter([0 -a(2:end)],1,segment_sig);    % Estimated signal
    e = segment_sig - est_x;
    LPE = xcorr(e,'coeff');   % LPES
    LPE=LPE(length(LPE)/2:end);

    for ii=1:length(hd)
        sigBW(:,ii)=filter(hd(ii),LPE);
        sig_TKEO(ii) = mean(TKEO(sigBW(:,ii)));
        sig_energy(ii) = mean(sigBW(:,ii)).^2; 
    end
    Hilb_tr = hilbert(sigBW);
    Hilb_env = abs(Hilb_tr);
    c = xcorr(Hilb_env);
    [cval,cidx] = max(c);
    GNEm(i) = max(cval);
    
    signal_BW_TKEO(i,:) = sig_TKEO;  
    signal_BW_energy(i,:) = sig_energy;
end

signal_BW_TKEO2 = mean(log(signal_BW_TKEO)); % used for getting the noise to signal ratio
signal_energy2 = mean(log(signal_BW_energy));

% Set outputs

GNE(1) = mean(GNEm);
GNE(2) = std(GNEm);

gnTKEO = mean(signal_BW_TKEO);
gnSEO = mean(signal_BW_energy);
GNE(3) = sum(gnTKEO(1:2))/sum(gnTKEO(end-3:end));
GNE(4) = sum(gnSEO(1:2))/sum(gnSEO(end-3:end));
GNE(5) = sum(signal_BW_TKEO2(end-3:end))/sum(signal_BW_TKEO2(1:2));
GNE(6) = sum(signal_energy2(end-3:end))/sum(signal_energy2(1:2));

end

function [VFER] = VFER_measure(data, fs, VF_close, VF_open)

filt_order=100;
BW = 500; % bandwidth 
Fmax = (fs/2-BW-300); % Max frequency to check 
Fshift = 250; % shift 

% cut-off1
Fc1=1:Fshift:Fmax;
% cut-off2
Fc2=Fc1+BW;

for j=1:length(Fc1);
    d(j) = fdesign.bandpass('n,fc1,fc2',filt_order,Fc1(j),Fc2(j),fs);
    hd(j) = design(d(j));
end

for i=1:length(VF_close)-1
    tseries = data(VF_close(i):VF_close(i+1));
    Dwindow = hann(length(tseries)); %Use Hanning window
    segment_sig = tseries.*Dwindow;
    
    if (length(tseries)>50)
        for ii=1:length(hd)
            thanasis = filter(hd(ii),segment_sig);
            sigBW(:,ii) = thanasis(1:50);
            sig_TKEO(ii) = mean(TKEO(sigBW(:,ii)));
            sig_SEO(ii) = mean(sigBW(:,ii)).^2; 
        end
        Hilb_tr = hilbert(sigBW);
        Hilb_env = abs(Hilb_tr);
        c = xcorr(Hilb_env);
        [cval,cidx] = max(c);
        NEm(i) = max(cval);

        signal_BW_TKEO(i,:) = sig_TKEO;
        signal_BW_SEO(i,:) = sig_SEO;        
    end
end

signal_BW_TKEO2 = mean(log(signal_BW_TKEO)); % used for getting the noise to signal ratio

% Set outputs

VFER(1) = mean(NEm);
VFER(2) = std(NEm);
VFER(3) = -sum(NEm.*log_bb(NEm));
VFTKEO = mean(signal_BW_TKEO);
VFSEO = mean(signal_BW_SEO);
VFlog_SEO = mean(log(signal_BW_SEO));

% Get 'signal to noise' ratios
size(signal_BW_TKEO)
VFER(4) = sum(VFTKEO(1:5))/sum(VFTKEO(6:10));
VFER(5) = sum(VFSEO(1:5))/sum(VFSEO(6:10));
VFER(6) = sum(signal_BW_TKEO2(6:10))/sum(signal_BW_TKEO2(1:5));
VFER(7) = sum(VFlog_SEO(6:10))/sum(VFlog_SEO(1:5));

end



function [HNR, NHR] = HNRFun(data,fs)

f0max=500; %Hz -- max value, possibly adjust for other applications
f0min=50; %Hz
tstep=0.01*fs;
x=0.08*fs;
steps=(length(data)-x)/tstep;

for i=1:steps
    
    tseries = data(i*tstep:i*tstep+x);
    tseries = tseries-mean(tseries);
    Dwindow = hann(length(tseries));
    segment_sig = tseries.*Dwindow;
    
    %% HNR computation process
    ACF = xcorr(segment_sig,'coeff');
    ACF2 = ACF(length(segment_sig):end);
    aa=fft(segment_sig);
    aa=ifft(abs(aa).^2);
    ACF_Dwindow = xcorr(Dwindow,'coeff');
    ACF_Dwindow2 = ACF_Dwindow(length(Dwindow):end);
    bb=fft(Dwindow);
    bb=ifft(abs(bb).^2);
    ACF_signal = ACF2./ACF_Dwindow2;
    ACF_signal = ACF_signal(1:round(length(ACF_signal)/3));
    rho=aa./bb;
    rho=rho(1:length(rho)/2);
    rho=rho/max(rho);
    [rx_value,rx_index] = sort(ACF_signal,'descend');
    [d1 d2] = sort(rho, 'descend');
    low_lim=ceil(fs/f0max);  % round towards positive sample number
    up_lim=floor(fs/f0min);  % round towards negative sample number
    k=2;
    while ((rx_index(k)<low_lim) || rx_index(k)>up_lim)
        k=k+1;
    end
    
    m=2;
    while ((d2(m)<low_lim) || d2(m)>up_lim)
        m=m+1;
    end
    ll(i)=d2(m);
    mm=d2(m); 
    HNR_dB_Praat(i) = 10*log10(rho(mm)/(1-rho(mm)));
    NHR_Praat(i) = (1-rho(mm))/rho(mm);

end

%% Summarize data
HNR(1)=mean(HNR_dB_Praat);
HNR(2)=std(HNR_dB_Praat);

NHR(1)=mean(NHR_Praat);
NHR(2)=std(NHR_Praat);

end



function pout = log_bb(pin, method)
% Function that computes the algorithm depending on the user specified
% base; if the input probability is zero it returns zero.

if nargin<2
    method = 'Nats';
end

switch (method)
    case 'Hartson' % using log10 for the entropy computation
        log_b=@log10;
        
    case 'Nats' % using ln (natural log) for the entropy computation 
        log_b=@log;
       
    otherwise % method -> 'Bits' using log2 for the entropy computation 
        log_b=@log2;
end

if pin==0
    pout=0;
else
    pout=log_b(pin);
end

end


function H = H_entropy(f)
% Calculate entropy of a discrete distribution
% Usage: H = entropy(f)
%  f - input distribution as a vector
%  H - entropy
N = length(f);
H = 0;
for j = 1:N
   H = H - f(j) * log_bb(f(j));
end

end % end of H_entropy function
