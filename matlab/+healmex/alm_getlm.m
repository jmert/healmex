function [l,m] = alm_getlm(lmax,idx)
% [l,m] = alm_getlm(lmax,idx)
%
% Get the l and m from index and lmax.
% lmax : The maximum l associated to the alm
% idx : The index for which to compute the l and m.

  if ~exist('idx', 'var')
    idx=1:healmex.alm_getn(lmax);
  end
  
  m = round(ceil(((2 * lmax + 1) - sqrt((2 * lmax + 1) ^ 2 - 8 * (idx - 1 - lmax))) / 2));
  l = idx - 1 - round(m .* (2 * lmax + 1 - m) / 2);
return