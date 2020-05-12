
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