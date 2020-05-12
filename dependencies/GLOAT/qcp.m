function [ LSF_vt, LSF_g, dglot,dglot2, a_vt] = qcp( frame,p_vt,p_g, gci_ins, fs )
%QCP Summary of this function goes here
%   Detailed explanation goes here

f2 = filter([1 -1],1,frame);



    
    PQ = 0.005;
    DQ = 0.7;
    d = 0.00001;
    %Nramp = round(0.000875*fs);
    %Nramp = round(0.0005*fs);
    Nramp = 1;
    w = makeW(frame,p_vt,DQ,PQ,d,Nramp,gci_ins,fs);
    wh = hann(length(f2));
    wh2 = hann(length(f2)-p_vt);
    
    a_vt = wlp(f2.*wh,w,p_vt);
    
        
if 0
    dglot2 = filter(a_vt,1,frame);
    a_g1 = lpc(dglot2(p_vt+1:end).*wh2,4);
    f3 = filter(a_g1,1,frame);
    a_vt = wlp(f3.*wh,w,p_vt);
    
end
    
    
    %a_vt = fixzeroes(a_vt,0,fs);
  
    
    if ~isminphase(a_vt)
      %  a_vt = swlp(f2,w,p_vt);
       % disp('HEP');
    end
    %a_vt2 = swlp(f2,w,p_vt);
    dglot = filter(a_vt,1,frame);
    dglot2 = dglot;
    %dglot2 = filter(a_vt2,1,frame);
    % figure(1);
    % clf;
    % zplane(a_vt,1);
    % waitforbuttonpress;

  if ~isminphase(a_vt)
        a_vt = a2minphase(a_vt,512);
  end

dglot = dglot(p_vt+1:end);
%dglot2 = dglot2(p_vt+1:end);
alpha = 0.99999;

glot = filter(1,[1 -alpha],dglot);

%a_g = lpc(glot,p_g);
a_g = lpc(hann(length(dglot)).*dglot,p_g);


LSF_vt = poly2lsf(a_vt);


%LSF_vt = poly2rc(a_vt);
%LSF_vt = rc2lar(LSF_vt);
LSF_g = poly2lsf(a_g);

gci_ins = gci_ins-p_vt;
gci_ins = gci_ins(gci_ins > 0);

%
% figure(1);
% clf
% hold on;
% plot(frame(p_vt+1:end)./norm(frame(p_vt+1:end),inf));
% plot(gci_ins,dglot(gci_ins),'xg','LineWidth',2);
% plot(dglot/norm(dglot,inf),'r');
% plot(w(p_vt+1:end),'m');
% waitforbuttonpress;



end

