function [ f0vec, mbf0 ] = estimate_f0(data,fs,wlen,Nhop,method)
%ESTIMATE_F0 Estimates the fundamental frequency of the given LPC residual
%signal by the autocorrelation method. Returns 0 if the frame is determined
%unvoiced.
%   
   Ndata = length(data);
    Nframes = floor((Ndata-wlen)/Nhop);

f0vec = zeros(Nframes,1);

%disp('Estimating F0...');

% COMPUTE PARAMETERS


for i = 1:Nframes

  start = 1+(i-1)*Nhop;
  stop = start+wlen-1;
  frame = data(start:stop);
  if strcmp(method,'yin')
    f0vec(i) = yin(frame,fs);
  else
      f0vec(i) = ac(frame,fs);
  end
  
end


Nw = round(0.5*fs/Nhop);

mbf0 = mbsig(f0vec,Nw);

end



function y = mbsig(x,Nw)
x2 = [zeros(Nw,1); x; zeros(Nw+1,1)];
%w = blackman(2*Nw+1);
w = ones(2*Nw+1,1);
N = length(x);
y = zeros(N,1);

for n = 1:N
    for m=-Nw:Nw
        y(n) = y(n) + w(m+Nw+1)*x2(n+m+Nw+1);    
    end
        %y(n) = y(n)/(2*Nw+1);
end
    y = y/sum(w);
end

function f0 = yin(frame,Fs)
    %frame = filter(lpc(frame,18),1,frame);
    W = length(frame);
    frame = frame(:);
    %%%%%%%%%%%% YIN Algorithm
    
    
%%%% Step 1: the autocorrelation funtcion:  
    dt = zeros(W,1);
    frame = [frame;zeros(W,1)];
    
    for tau = 1:W
        for j = 1:W
            dt(tau) = dt(tau) + (frame(j)-frame(j+tau)).^2;
        end
    end
    
    %%%% Step 3: Cumulative mean normalized difference function
    
    dt2 = zeros(W,1);
    dt2(1) = 1;
    for tau = 2:W
        summa = 0;
        for j = 2:tau
            summa = summa + dt(j);
        end
        dt2(tau) = dt(tau)/(summa/(tau-1));
    end;
    
    %%%% Step 4: Absolute threshold
    minIndex = round(Fs/500);
    maxIndex = round(Fs/50);
    dt3 = dt2(minIndex:maxIndex);
    threshold = 0.1;
    I = find(dt3 <= threshold);
    if isempty(I)
        [Y ind] = min(dt3);
    else
        ind = I(1);
        range = ind:ind+5;
        [Y ind] = min(dt3(range));
    end
    
    %%%% Step 5: Parabolic interpolation
    
    if ind > 1 && ind < length(dt3);
    alpha = dt3(ind-1);
    beta = dt3(ind);
    gamma = dt3(ind+1);
    p = 0.5*(alpha-gamma)/(alpha-2*beta+gamma);
    else
        p = 0;
    end
    
    f0 = Fs/(ind+minIndex-1+p-1);
   
end  
   
function f0 = ac(frame,Fs)    
 %%%%%%%%%%%% Autocorrelation Method    
    p = round(Fs/1000+2);
    frame = frame(:);
    w = hann(length(frame));
    
    frame_lpc = lpc(w.*frame, p);
    lpc_residual = filter(frame_lpc,1,frame);
    
    Fr = fft(lpc_residual,length(lpc_residual)*2);
    Sf = Fr.*conj(Fr);
    R = real(ifft(Sf)); % Autocorrelation function
    
    tau_max = round(Fs/50);
    tau_min = round(Fs/350);
    
    r_max = 0;
    
    for i = tau_min:tau_max
        if abs(R(i)) > r_max
            r_max = abs(R(i));
            f0 = Fs/(i+1);
        end
    end
    


end

