function l = len(x)
% LEN  Length of the signal

% $Id: len.m 3 2004-02-04 12:57:04Z mairas $

for i=1:length(x)
  l(i) = length(x(i).s);
end
