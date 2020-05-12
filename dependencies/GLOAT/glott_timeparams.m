function q=glott_timeparams(g,f0,qoq_level)
% Taken from Aparat and slightly modified

% GLOTTALTIMEPARAMS Time-based parameterisation of the glottal flow
% Q=GLOTTALTIMEPARAMS(G,F0) returns a structure containing
% multiple time-based parameters of the glottal flow together with
% the time instants required for their calculation.
%
% Description of the different fields:
%
% OQ1
% OQ2
%
% OQ1 and OQ2 are the open quotients calculated from the so-called primary
% and secondary openings of the glottal flow, respectively. Please see
% Hannu Pulakka's master's thesis
% (http://www.acoustics.hut.fi/publications/files/theses/pulakka_mst.pdf)
% around pages 37-38 to see a description of the two openings. If unsure,
% use OQ1.
%
% NAQ
%
% The Normalized Amplitude Quotient, a comparably robust parameter closely
% related to the pressedness of speech. For more information, see:
%
% P. Alku, T. Bäckström, and E. Vilkman. Normalized Amplitude Quotient for
% Parameterization of the Glottal Flow. Journal of the Acoustical Society
% of America, 112(2):701-710, August 2002.
%
% AQ
%
% The Amplitude Quotient:
%
% P. Alku and E. Vilkman. Amplitude Domain Quotient for Characterization
% of the Glottal Volume Velocity Waveform Estimated by Inverse Filtering.
% Speech Communication, 18(2):131-138, 1996.
%
% ClQ: 0.1751
%
% The Closing Quotient.
%
% OQa
%
% A variation of the Open Quotient derived from the LF model:
%
% Gobl, C. and Ní Chasaide, A., 2003. Amplitude-based source parameters
% for measuring voice quality Proceedings of the ISCA Tutorial and
% Research Workshop VOQUAL'03 on Voice Quality: Functions, Analysis and
% Synthesis, Geneva, pp. 151-156.
%
% QOQ
%
% The Quasi-Open Quotient. See e.g.:
%
% A.-M. Laukkanen, E. Vilkman, P. Alku, and H. Oksanen. Physical
% Variations Related to Stress and Emotional State: a Preliminary Study.
% Journal of Phonetics, 24(3):313-335, 1996.
%
% SQ1
% SQ2
%
% Speed Quotients calculated from the primary and the secondary openings,
% respectively (see OQ1 & OQ2).
%
% t_po
%
% The time instant of the primary opening.
%
% t_so
%
% The time instant of the secondary opening.
%
% t_c
%
% The closing time instant.
%
% t_max
%
% The time instant of the maximum flow.
%
% t_min
%
% The time instant of the minimum flow.
%
% t_dmin
%
% The time instant of the minimum of the derivative.
%
% t_dmax
%
% The time instant of the maximum of the derivative.
%
% t_qo
%
% The quasi-opening time instant.
%
% t_qc
%
% The quasi-closing time instant.
%

% $Id: glottaltimeparams.m 200 2007-09-10 12:32:40Z mairas $

if nargin<2
  f0 = find_f0_yin(g)
end
if nargin<3
  qoq_level = 0.5;
end

% set the primary opening detection threshold
lvl = 0.2;

n0 = round(g.fs/f0);
T0 = 1/f0;

% smoothen the waveform

b_sm = fir1(4,4000/g.fs/2);

%g_sm = filter_ap(b_sm,1,g,'noncausal'); %% commented
b = mk_gauss_filter(T0*g.fs);
g_mex = filter(b,1,g.s); %,'noncausal'); %% commented 

g_sm = g;
% Get the absolute maximum of the frame. It is supposed to be one
% of the Amax values.

mx1 = max(g_sm.s);
[pks, I_gmax] = findpeaks(g_sm.s,'MinPeakDistance',n0*0.4,'MinPeakHeight',mx1*0.50);

%[gmax,I_gmax]=max(g_sm.s); 
t_gmax=g_sm.t(I_gmax);

% Get the time values of the maxima
%T=repeatmax(t_gmax,f0,g_sm);
T = I_gmax; %t_gmax;

% init q
% q=struct('OQ1', {}, 'OQ2', {}, ...
%   'NAQ', {}, 'AQ', {}, ...
%   'ClQ', {}, 'OQa', {}, ...
%   'QOQ', {}, ...
%   'SQ1', {}, 'SQ2', {}, ...
%   't_po', {}, 't_so', {}, ...
%   't_c', {}, 't_max', {}', ...
%   't_min', {}, 't_dmin', {}, ...
%   't_dmax', {}, 't_qo', {}, 't_qc', {});

q = [];
j=1;
for i=1:length(T)-1
    
    t_max=T(i);
    %T
   %A_max=g_sm.at(t_max);
   A_max=g_sm.s(round(t_max));
   
   [g_per,g_mex_per] = get_period_ap(g_sm,t_max,n0,g_mex);
%    n0
%    length(g_per)
   
   if length(g_per) == 1
     if g_per==1
       break;
     elseif g_per==-1
       continue;
     end
   end
  
  
  dg_per = diff(g_per); dg_per(length(dg_per)+1) = dg_per(end); %% deriv(g_per);
  
  [~,gper_max] = max(g_per);
  if gper_max > n0
      continue;
  end
  
  g_p = g_per(1:gper_max);
  g_n = g_per(gper_max:end);
  
  dg_p = dg_per(1:gper_max);
  dg_n = dg_per(gper_max:end);
  
  [d_max,I_tdmax]=max(dg_p(round(0.05*n0):end));
  
  [A_min,I_tmin]=min(g_n);
  
  [d_min,I_tdmin]=min(dg_n);
  d_min=-d_min; % get the amplitude of the peak
  
  %t_min=g_n.t(I_tmin);
  t_min = I_tmin + gper_max -1;
  t_dmin = I_tdmin + gper_max -1;
  
  t_dmax = round(0.05*n0) + I_tdmax -1;
  
  A_ac=A_max-A_min;

  l = A_min+A_ac*lvl;
  
  % get the positive threshold level crossing point
  
  %t_p=crossings(g_p,l,1);
  idx=find(diff((g_p - l)>0)>0)+1;
%   if (isempty(idx)) %% check this 
%     idx = 1;
%   end
  t_p = idx;
  
  
  if length(t_p)==0
      continue
  end
  t_p=t_p(end);
  
  % get the first local minimum before the 10% amplitude threshold
  g_s = g_per(1:t_p(1));
  
  idx = find(diff(( diff(g_s)  -0)>0)>0)+1;
  t_cr = idx;
    
  if length(t_cr)==0
    t_po=1; % g_s.t(1);
  else
    t_po = t_cr(end);
  end
  A_po = g_s(t_po);
  
  % get the location of the primary opening
  t_po1 = t_po;
  while 1
      if ((t_po-round(0.05*n0)) > 1) 
          g_s = g_per(t_po-round(0.05*n0):t_po); %trim(g_per,t_po-0.05*T0,t_po);
          t_po1 = t_po-round(0.05*n0);
      else
          g_s = g_per(t_po:t_po);
          t_po1 = t_po;
      end
    [y,I] = min(g_s);
    if y<A_po-0.01*A_ac
      t_po = I + t_po1 -1; % g_s.t(I);
      A_po = g_s(I);
      break;
    else
      break;
    end
  end
%   t_po1
%   t_po
  % get the location of the secondary opening
  %g_s = g_per(round(t_po+0.05*n0):gper_max); %trim(g_mex,t_po+0.05*T0,t_max);
  
  if round(t_po+0.05*n0) < gper_max
      g_s = g_mex_per(round(t_po+0.05*n0):gper_max);
  else
      g_s = g_mex_per(round(t_po+0.05*n0):round(t_po+0.05*n0));
  end
  [c,I]=max(g_s);
  t_so = round(t_po+0.05*n0) + I -1;    %g_s.t(I);
  
  % get the first positive derivative zero-crossing after the first minimum
  % slope after t_max
  
  g_n2 = g_n(I_tdmin:length(g_n));
  idx = find(diff(( diff(g_n2) - 0)>0)>0)+1;
  
  t_dzc = gper_max + I_tdmin + idx -2; %crossings(diff(g_n2),0,1);
  
  if length(t_dzc)==0
    % if no local minima are found, we must be at the end of the
    % signal.
    break;
  end

  % get the QOQ positive threshold level crossing poing
  qoql = A_min+A_ac*qoq_level;
  t_ps = find(diff(( g_p  - qoql)>0)>0)+1; %  crossings(g_p,qoql,1);
  if length(t_ps)==0, continue, end
  t_qo=t_ps(end);
  % get the positive threshold level crossing point

  t_ns = find(diff(( g_n  - qoql)>=0)<0)+1;  % crossings(g_n,qoql,-1);
  if length(t_ns)==0, continue, end
  t_qc= gper_max + t_ns(1)-1;
  
  t_c = t_dzc(1);
  T_c = t_c-gper_max; % t_max;
  T_po = t_c-t_po;
  T_so = t_c-t_so;
  T_pop = gper_max-t_po;  % t_max-t_po;
  T_cl = t_c-gper_max;    % t_max;
  T_sop = gper_max-t_so;  %t_max-t_so;
  T1 = t_qc-t_qo;
  
  ClQ = T_c/n0; % T0;
  AQ = A_ac/d_min;
  NAQ = AQ/n0;  % T0;
  
  OQ1 = T_po/n0; %T0;
  OQ2 = T_so/n0; %T0;
  
  OQa = A_ac*(pi/(2*d_max)+1/d_min)*f0;
  
  QOQ = T1/n0; %T0;
    
  SQ1=T_pop/T_cl;
  SQ2=T_sop/T_cl;
  
  t_max = gper_max; %% added to make t_max relative
  % store the values in a structure array
%   q(j)=struct('OQ1', OQ1, 'OQ2', OQ2, ...
% 	      'NAQ', NAQ, 'AQ', AQ, ...
% 	      'ClQ', ClQ, 'OQa', OQa, ...
% 	      'QOQ', QOQ, ...
% 	      'SQ1', SQ1, 'SQ2', SQ2, ...
% 	      't_po', t_po, 't_so', t_so, ...
% 	      't_c', t_c, 't_max', t_max', ...
% 	      't_min', t_min, 't_dmin', t_dmin, ...
% 	      't_dmax', t_dmax, 't_qo', t_qo, 't_qc', t_qc);
     q(j,1) = OQ1;
     q(j,2) = OQ2;
     q(j,3) = NAQ;
     q(j,4) = AQ;
     q(j,5) = ClQ;
     q(j,6) = OQa;
     q(j,7) = QOQ;
     q(j,8) = SQ1;
     q(j,9) = SQ2;
     q(j,10) = t_c;
        
  j=j+1;
end

1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=repeatmax(t0,f0,x)

% number of periods fitting the region
N = floor((t0-x.t(1))*f0)+floor((x.t(end)-t0)*f0)+1;

% index of the given location
n0 = floor((t0-x.t(1))*f0)+1;

T=[];

for i=1:N
  T = [T maxnear(x, t0+(i-n0)/f0, 1/f0/2)];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t=maxnear(x,t0,r)
% MAXNEAR Find the position of local maximum of X near T0 (search radius R)

x_ = trim_ap(x,t0-r,t0+r);

[x_m,I] = max(x_);
t = x_.t(I);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = mk_gauss_filter(L0)

% create a mexican hat filter. L0 is the length of a period, in samples.

gaussian_r = round(L0/20);
win_r = 5*gaussian_r;
t = -win_r:win_r;
b = exp(-(t./gaussian_r).^2);
b = diff(b,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [g_per,g_mex_per] = get_period_ap(g,t_max,n0,g_mex)

% acquire a single period of the signal

% returns either:
%   the period as a signal object
%   -1 if at the beginning of the signal
%   +1 if at the end of the signal

% get the next minimum slope

if abs(length(g.s)-t_max)<1e-10
  g_per = 1;
  g_mex_per = 1;
  return;
end
if (t_max+n0) > length(g.s)
  g_per = 1;
  g_mex_per = 1;
  return;
end


g_next = g.s(round(t_max):round((t_max+n0)));
d_gnext=diff(g_next);
[d_min,idmin]=min(d_gnext);
t_dmin=(idmin);

% get the first positive derivative zero-crossing after that
  
g_next2=g.s(t_dmin+t_max:t_max+n0);
tm = find(diff(( diff(g_next2)  - 0)>0)>0)+1;  %find(diff(g_next2)>0);
%tm=crossings(diff(g_next2),0,1);


if length(tm)==0
  % if no local minima are found, we must be at the end of the
  % signal.
  g_per = 1;
  g_mex_per = 1;
  return;
end
  
t_c = tm(1)+t_dmin+t_max-1;
  
% define pulse begin and end times
  
t_pb = t_c-n0;
t_pe = t_c+round(0.2*n0);
  
if t_pb < 1 %g.t(1)
  % need more margin at the beginning of the signal
  g_per = -1;
  g_mex_per = -1;
  return;
end
  
if t_pe > length(g.s)  %g.t(end)
  % need more margin at the end of the signal
  g_per = 1;
  g_mex_per = 1;
  return;
end
  
%g_per = trim(g,t_pb,t_pe);
g_per = g.s(t_pb:t_pe);
g_mex_per = g_mex(t_pb:t_pe);


function y = deriv(x,n,l,varargin)
% DERIV(X) Calculate the approximate derivative of the signal.
% DERIV(X,N) n'th order approximation (default N=1)
% DERIV(X,N,L) l'th derivative (default L=1)

% $Id: deriv.m 120 2007-01-02 12:21:49Z mairas $

if nargin < 3
    if nargin < 2
        n = 1;
    end
    l = 1;
elseif length(n) == 0
    n = 1;
end

k = n/2;
s = dsinc(-k:k).*hanning(k*2+1)';
s = -s./(s*(0:n)');

y = conv(s,x.s);
y = y(1+ceil(n):floor(length(y)-n));
t = x.time;
t = set(t,'num', length(y), 'begin', t.begin+(k-floor(k))/t.fs);
y = signal(y,t);
y = y*y.time.fs;

if l > 1 % derivative recursion - not neat but it works
    y = deriv(y,n,l-1);
end


% Old version 
%y = diff(x,varargin{:});
%y = y*y.time.fs;



function y = dsinc(x)
% function y = dsinc(x)
% The sinc function derivative

i=find(x==0);
x(i)= 1;      
y = cos(pi*x)./(pi*x) - sin(pi*x)./((pi*x).^2);
y(i) = 0;
 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g_per = get_period(g,t_max,T0);

% acquire a single period of the signal

% returns either:
%   the period as a signal object
%   -1 if at the beginning of the signal
%   +1 if at the end of the signal

% get the next minimum slope

if abs(g.time.t(end)-t_max)<1e-10
  g_per = 1;
  return;
end

g_next=trim_ap(g,t_max,t_max+T0);
d_gnext=diff(g_next);
[d_min,idmin]=min(d_gnext);
t_dmin=d_gnext.t(idmin);

% get the first positive derivative zero-crossing after that
  
g_next2=trim(g,t_dmin,t_max+T0);
tm=crossings(diff(g_next2),0,1);
  
if length(tm)==0
  % if no local minima are found, we must be at the end of the
  % signal.
  g_per = 1;
  return;
end
  
t_c = tm(1);
  
% define pulse begin and end times
  
t_pb = t_c-T0;
t_pe = t_c+0.2*T0;
  
if t_pb < g.t(1)
  % need more margin at the beginning of the signal
  g_per = -1;
  return;
end
  
if t_pe > g.t(end)
  % need more margin at the end of the signal
  g_per = 1;
  return;
end
  
g_per = trim(g,t_pb,t_pe);
