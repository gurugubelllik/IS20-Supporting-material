function [vnvsig,vnvevi,zf,gci,es,f0] = zff_analysis(s,fs,nmean,VNVZFTH)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage:  [vnvsig,vnvevi,zf,gci,es,f0]  =  zff_analysis(s,fs,nmean,VNVZFTH)
% Authors: Dhananjaya N
% Algorithm: Sri Rama Murty and Professor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modifed by:
%        G Krishna  (Speech Lab, IIIT Hyderabad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~exist('nmean'))
 nmean=10;
end;

s=s(:);
ds=diff(s);
ds(end+1)=ds(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r = 0.98;
% f = 0;
% [B1,A1]=resonator(r,f,fs);
% 
% zf=filter(B1,A1,ds);
% zf=filter(B1,A1,zf);
%
% musig=RunMean(zf,floor(nmean*fs/1000));
% zf=zf-musig;
% 
% musig=RunMean(zf,floor(nmean*fs/1000));
% zf=zf-musig;
% 
% musig=RunMean(zf,floor(nmean*fs/1000));
% zf=zf-musig;
% 
% musig=RunMean(zf,floor(nmean*fs/1000));
% zf=zf-musig;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[zf,gci,es,f0]=zfsig(ds,fs);
zfe=RunMean(zf.^2,nmean*fs/1000);
zfe=zfe/max(zfe);
zfevi=1-exp(-10*zfe);

if(~exist('VNVZFTH'))
    VNVZFTH=1-exp(-10*.005);
else
    VNVZFTH=1-exp(-10*VNVZFTH);
end

vnvevi=zfevi;
vnvsig=vnvevi > VNVZFTH;
vnvsig=RunMean(double(vnvsig),round(nmean*fs/1000));
vnvsig=vnvsig>VNVZFTH;

% figure;clear ax;
% ax(1)=subplot(3,1,1);plot(s);
% ax(2)=subplot(3,1,2);plot(zf);
% ax(3)=subplot(3,1,3);plot(gci,es/max(es),'.');hold on;plot(zfe,'r');
% linkaxes(ax,'x');

return;

end




function [zf,gci,es,f0,extrainfo] = zfsig(wav,fs,winLength,wintype)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [zf,gci,es,f0] = zfsig(wav,fs,winLength)
%
%   Computes the zero-frequency filtered signal for the input speech signal
%   'wav' at sampling rate 'fs'. 'winLength' is an optional parameter in ms
%   which specifies the window length for trend removal.
%
% Returns:
%       zf      - zero-frequency resonator signal
%       gci     - glottal closure instants
%       es      - excitation strength at GCIs
%       f0      - pitch frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~exist('winLength') | winLength == 0)
    winLength=xcorrWinLen(wav,fs);
    %winLength=1.2*winLength;
end

if(~exist('wintype'))
    wintype='rec';
end

zf=zeroFreqFilter(wav,fs,winLength,wintype);
%figure;plot(wav);pause;

pgci=find(diff((zf>0))==1); % +ve zero crossings
ngci=find(diff((zf>0))==-1); % +ve zero crossings

pes=abs(zf(pgci+1)-zf(pgci-1));
nes=abs(zf(ngci+1)-zf(ngci-1));

avgpes=mean(pes);
avgnes=mean(nes);

if(avgpes>avgnes)
    gci=pgci;
    es=pes;
else
    gci=ngci;
    es=nes;
end
T0=diff(gci);
T0=T0(:)/fs;
f0=1./T0;
%f0(end+1)=f0(end);
f0=[f0(1);f0];

extrainfo.pgci=pgci;
extrainfo.ngci=ngci;
extrainfo.pes=pes;
extrainfo.nes=nes;
extrainfo.winlen=winLength;
end

function [zfSig]=zeroFreqFilter(wav,fs,winLength,wintype)

dwav=diff(wav);
dwav=[dwav(1);dwav(:)];
dwav=dwav/max(abs(dwav));
N=length(dwav);

%zfSig=filter(1,[1 -4 6 -4 1],dwav);
zfSig=cumsum(cumsum(cumsum(cumsum(dwav))));

winLength=round(winLength*fs/1000);
zfSig=remTrend(zfSig,winLength,wintype);
zfSig=remTrend(zfSig,winLength,wintype);
zfSig=remTrend(zfSig,winLength,wintype);
%zfSig=remTrend(zfSig,winLength,wintype);
zfSig(N-winLength*2:N)=0;
zfSig(1:winLength*2)=0;
end

function [out]=remTrend(sig,winSize,wintype)

switch(lower(wintype))
    case {'rec','rect','rectangle','1'}
        window = ones(winSize,1);
    case {'kai','kaiser','2'}
        a = 50;
        b = 0.1102*(a-8.7);
        window = kaiser(winSize,b);
        %window=window/sqrt(sum(window.^2));
        window = window/sum(window)*winSize;
    case {'HAM', 'HAMM', 'ham', 'hamm', '2'}
        window = hamming(winSize);
        window = window/sum(window)*winSize;
    case {'HAN', 'HANN', 'han', 'hann', '3'}
        window = hanning(winSize);
        window = window/sum(window)*winSize;
    otherwise
        error('Wrong wintype...');
end
%window(winSize/2-3:winSize/2+3)=0;

rm=conv(sig,window);
rm=rm(ceil(winSize/2):length(rm)-(winSize-ceil(winSize/2)));

norm=conv(ones(size(sig)),window);
norm=norm(ceil(winSize/2):length(norm)-(winSize-ceil(winSize/2)));

rm=rm./norm;
out=sig-rm;
end

function [idx]=xcorrWinLen(wav,fs)

% 	zfSig=zeroFreqFilter(wav,fs,2);
% 	zfSig=zfSig/max(abs(zfSig));
% 	wav=zfSig;

frameSize=30*fs/1000;
frameShift=20*fs/1000;

en=conv(wav.^2,ones(frameSize,1));
en=en(frameSize/2:end-frameSize/2);
en=en/frameSize;
en=sqrt(en);
en=en>max(en)/5;

b=buffer(wav,frameSize,frameShift,'nodelay');
vad=sum(buffer(en,frameSize,frameShift,'nodelay'));

FUN=@(x) xcorr((x-mean(x)).*hamming(length(x)),'coeff')./xcorr(hamming(length(x)),'coeff');
out=blkproc(b,[frameSize,1],FUN);

out=out(frameSize:end,:);

minPitch=3;  %2 ms == 500 Hz.
maxPitch=16; %16 ms == 66.66 Hz.

[maxv maxi]=max(out(minPitch*fs/1000:maxPitch*fs/1000,:));

%h=hist(maxi(vad>frameSize/2)+minPitch,(3:15)*8-4);
x=(minPitch:0.5:maxPitch)*fs/1000+2;
pLoc=maxi(vad>frameSize*0.8)+minPitch*fs/1000;
y=hist(pLoc,x);
y=y/length(pLoc);

%h=figure();
%bar(x,y,1,'EdgeColor',[1 1 1],'FaceColor',[0 0 0]);
%set(gca,'xTick',(1:maxPitch)*fs/1000+0.5*fs/1000, 'xTickLabel',(1:maxPitch));
%set(gca,'yTick',[0 0.1 0.2 0.3 0.4],'yTickLabel',[0 0.1 0.2 0.3 0.4]);
%xlabel('Time (s)');
%ylabel('Normalized frequency');
%allText   = findall(h, 'type', 'text');
%allAxes   = findall(h, 'type', 'axes');
%allFont   = [allText; allAxes];
%xlim([1 maxPitch+1]*fs/1000)
%set(allFont,'FontSize',18);

%advexpfig(h,'hist.eps','-deps2c','w',20,'h',20);

%close(h);

[val idx]=max(y);
idx=round(idx/2)+minPitch+2;
%disp(['Average pitch period: ' num2str(idx) ' ms']);
end

function [yavg, ysum] = RunMean(sig, N, wintype)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Usage : [yavg,ysum] = RunMean(sig, N, wintype)
%
% OUTPUTS	:
%	yavg	- running mean of the signal
%	ysum	- running sum of the signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if( ~ exist('wintype'))
    wintype	= 'RECT';
end

M	= length(sig);
Nby2	= floor(N/2);
Ntail	= N - Nby2 - 1;

switch lower(wintype)
    case { 'rec', 'rect', 'rectangle','rectangular', '1'}
        h	= ones(N,1);
        
    case {'ham', 'hamm', 'hamming' '2'}
        h	= hamming(N);
        
    case {'han', 'hann', 'hanning', '3'}
        h	= hanning(N);
        
    case {'kai', 'kaiser', '4'}
        a = 50;
        b = 0.1102*(a-8.7);
        h = kaiser(N,b);
        h = h/sqrt(sum(h.^2));
        
    case {'blk', 'black', 'blackman', '5'}
        h	= hanning(N);
        
    otherwise
        disp('Error : Unknown window type!!!');
        exit(1);
end

x	= conv(double(sig),double(h));
ysum	= x(Nby2+1:M+N-1-Ntail);

xdiv	= conv(ones(size(sig)),h);
x	= x ./ xdiv;
yavg	= x(Nby2+1:M+N-1-Ntail);

end