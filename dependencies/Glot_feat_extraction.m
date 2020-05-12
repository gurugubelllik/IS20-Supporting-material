function Glot_feat=Glot_feat_extraction(data,fs,vnv_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code to find Glott features from the speech Signal
% Authors:  G Krishna  (Speech Lab, IIIT Hyderabad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------       inputs         --------
% data      : input raw speech signal
% fs        : sampling freqnecy
% --------      outputs         --------
% I_feat    : voice quality features (192 dimentional)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dependencies
%       Functions in ./dependencies/GLOTT folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    vnv_flag=1;
end 

if(vnv_flag==1)
    [vnvsig,vnvevi,zf,gci_locs,es,f0 ]= zff_analysis(data,fs);    
    data=data(vnvsig==1);
end

wlen = 0.03*fs; % Frame length (sustained vowels, so long frame length is OK)
Nhop = 0.010*fs; % Hop size

if(fs==16000)
    p_vt = 24; % Vocal Tract filter order %24 for 16 Khz % 12 for 8 KHz
    p_vt2 = 24;
    p_g = 10; % Glottal Source model order %10 for 16khz  % 6 for 8 KHz
else
    p_vt = 12; % Vocal Tract filter order %24 for 16 Khz % 12 for 8 KHz
    p_vt2 = 12;
    p_g = 6; % Glottal Source model order %10 for 16khz  % 6 for 8 KHz
end

E_REF = 1; 

%data = filtfilt(hp,1,data);
Ndata = length(data);
Nframes = floor((Ndata-wlen)/Nhop);

%VUVref = vuv(data,wlen,Nhop,0.0,wName);
%VUVdata = interp1(linspace(0,1,Nframes),VUVref,linspace(0,1,Ndata),'linear');
%VUVdata(VUVdata >= 0.5) = 1;
%VUVdata(VUVdata < 1) = 0;
VUVdata = ones(Nframes,1);

[f0vec2, ~] = estimate_f0(data,fs,wlen,Nhop,'ac');

meanf0 = median(f0vec2);
% F0min = 100; F0max = 400;
% [f0_vec,VUVDecisions,SRHVal] = SRH_PitchTracking(data,fs,F0min,F0max);
% 
% if meanf0 < 150
%     meanf02 = 120;
% else
%     meanf02 = 180;
% end

gcis = [];

lsf_general = zeros(Nframes,p_vt2);
lsf_vt = zeros(Nframes,p_vt);
lsf_g = zeros(Nframes,p_g);

Evec = zeros(Nframes,1);
dglots = zeros(Nframes,wlen);
ifexc = zeros(1,length(data));
f0vec = zeros(Nframes,1);

sampleFrames = ones(size(data));

prevStop = 0;

for i = 1:Nframes
    if i > 2
        start = 1+(i-1)*Nhop;
        stop = start+wlen-1;
        start = start-p_vt;
    else
        start = 1;
        stop = start+wlen+p_vt-1;
    end
    frame = data(start:stop);
    Evec(i) = 10*log10(sum(frame(p_vt+1:end).^2)/E_REF/wlen);
    if Evec(i) < -30
        VUVdata(i) = 0;
    end

    [gci_ins, ~] = gci(frame,meanf0,fs);
    
    gci_ins2 = gci_ins+start-1;
    gci_ins2 = gci_ins2(gci_ins2 > prevStop);
    prevStop = stop;
    gci_ins2 = gci_ins2(:)';
    gcis = [gcis, gci_ins2];
    
    [LSF_vt, LSF_g, dglot,~] = qcp(frame,p_vt,p_g,gci_ins,fs);
    
    lsf_vt(i,:) = LSF_vt;
    lsf_g(i,:) = LSF_g;
    frame2 = filter([1 -1],1,frame);
    a = lpc(frame2,p_vt2);
    lsf_a = poly2lsf(a);
    lsf_general(i,:) = lsf_a;
    %Evec(i) = norm(frame(p_vt+1:end)-mean(frame(p_vt+1:end)))/wlen;
    Evec(i) = 10*log10(sum(frame(p_vt+1:end).^2)/E_REF/wlen);
    Eframe = norm(frame(p_vt+1:end)-mean(frame(p_vt+1:end)));
    dglot = dglot/norm(dglot-mean(dglot))*Eframe;
    dglots(i,:) = dglot;
    
    
    ifexc(start+p_vt:stop) = ifexc(start+p_vt:stop)+dglots(i,:).*hann(wlen)';
    start2 = round(stop-wlen/2-Nhop/2);
    stop2 = start2+Nhop;
    sampleFrames(start2:stop2) = i;
end
%% glott feat extraction

clear glt_time;
clear glt_freq;
j=1;
for i = 1:Nframes
    if i > 2
        start = 1+(i-1)*Nhop;
        stop = start+wlen-1;
        start = start-p_vt;
    else
        start = 1;
        stop = start+wlen+p_vt-1;
    end
    
%     if VUVdata(i) == 0
%         continue;
%     end
    
    
    dglot_f = ifexc(start:stop);
    
    nf0 = round((1/meanf0)*fs);
    mx1 = max(-dglot_f);
    [pks, gci_ins] = findpeaks(-dglot_f,'MinPeakDistance',nf0*0.4,'MinPeakHeight',mx1*0.50);
    
    
    if (length(gci_ins) < 3) || (std(diff(gci_ins)) > 20) || (gci_ins(1) > 0.6*length(dglot_f)) || (gci_ins(end) < 0.4*length(dglot_f))
        continue;
    end
    
    
    glot_f = filter(1,[1, -1.00],dglot_f);
    
    winSize = round(1.0*nf0);
    glot_f1 = remTrend(glot_f,winSize,'rec');
    glot_f1 = remTrend(glot_f1,winSize,'rec'); %%increases NAQ hence commented
    %gTar1 = gTar1-min(gTar1);
    glot_f = glot_f1;
    
    
    g.fs = fs;
    
    g.s = glot_f;
    t=1:length(g.s);
    t = t*(1/fs);
    g.t = t;
    
    g.time.fs = fs;
    g.time.begin = 0;
    g.time.end = length(g.s)*(1/fs);
    
    g.time.t = t;
    
    %g.at(g.t) = g.s;
    f0 = find_f0_yin(g);
    
    if (f0 < 80) || (f0 > 500)
        continue;
    end
    
    q = glott_timeparams(g,f0);
    
    if isempty(q)
        continue;
    end
    
    glt_time(j,1) = sum(q(:,1))/length(q(:,1));  %sum(q(1:end).OQ1)/length(q);
    glt_time(j,2) = sum(q(:,2))/length(q(:,2));  %sum(q(1:end).OQ2)/length(q);
    glt_time(j,3) = sum(q(:,3))/length(q(:,3));  %sum(q(1:end).NAQ)/length(q);
    glt_time(j,4) = sum(q(:,4))/length(q(:,4));  %sum(q(1:end).AQ)/length(q);
    glt_time(j,5) = sum(q(:,5))/length(q(:,5));  %sum(q(1:end).ClQ)/length(q);
    glt_time(j,6) = sum(q(:,6))/length(q(:,6));  %sum(q(1:end).OQa)/length(q);
    glt_time(j,7) = sum(q(:,7))/length(q(:,7));  %sum(q(1:end).QOQ)/length(q);
    glt_time(j,8) = sum(q(:,8))/length(q(:,8));  %sum(q(1:end).SQ1)/length(q);
    glt_time(j,9) = sum(q(:,9))/length(q(:,9));  %sum(q(1:end).SQ2)/length(q);
    
    
    q_f = glott_freqparams(g,f0);
    
    glt_freq(j,1) = q_f.DH12;
    glt_freq(j,2) = sum(q_f.PSP)/length(q_f.PSP);
    glt_freq(j,3) = q_f.HRF;
    
    j=j+1;
    
    if j == 1
        continue;
    end
    
    glt_mn = mean(glt_time); % 9
    glt_std = std(glt_time);
    glt_med = median(glt_time);
    glt_min = min(glt_time);
    glt_max = max(glt_time);
    glt_range = range(glt_time);
    glt_skew = skewness(glt_time);
    glt_kurt = kurtosis(glt_time);
    
    diff_glt_time = diff(glt_time);
    dglt_mn = mean(diff_glt_time);
    dglt_std = std(diff_glt_time);
    dglt_med = median(diff_glt_time);
    dglt_min = min(diff_glt_time);
    dglt_max = max(diff_glt_time);
    dglt_range = range(diff_glt_time);
    dglt_skew = skewness(diff_glt_time);
    dglt_kurt = kurtosis(diff_glt_time);
    
    glt_fr_mn = mean(glt_freq);
    glt_fr_std = std(glt_freq);
    glt_fr_med = median(glt_freq);
    glt_fr_min = min(glt_freq);
    glt_fr_max = max(glt_freq);
    glt_fr_range = range(glt_freq);
    glt_fr_skew = skewness(glt_freq);
    glt_fr_kurt = kurtosis(glt_freq);
    
    diff_glt_freq = diff(glt_freq);
    dglt_fr_mn = mean(diff_glt_freq);
    dglt_fr_std = std(diff_glt_freq);
    dglt_fr_med = median(diff_glt_freq);
    dglt_fr_min = min(diff_glt_freq);
    dglt_fr_max = max(diff_glt_freq);
    dglt_fr_range = range(diff_glt_freq);
    dglt_fr_skew = skewness(diff_glt_freq);
    dglt_fr_kurt = kurtosis(diff_glt_freq);
    
    Glot_feat = [glt_mn glt_std glt_med glt_min glt_max glt_range glt_skew glt_kurt  dglt_mn dglt_std dglt_med dglt_min dglt_max dglt_range dglt_skew dglt_kurt               glt_fr_mn glt_fr_std glt_fr_med glt_fr_min glt_fr_max glt_fr_range glt_fr_skew glt_fr_kurt  dglt_fr_mn dglt_fr_std dglt_fr_med dglt_fr_min dglt_fr_max dglt_fr_range dglt_fr_skew dglt_fr_kurt ];
    
end

end