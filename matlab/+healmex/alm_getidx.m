function idx = alm_getidx(lmax, l, m)
% idx = alm_getidx(lmax, l, m)
%
% Returns the index of the alm coefficient (l,m) within a vector which includes
% coefficients up to lmax.

  if ~isscalar(lmax)
    error('lmax must be a scalar');
  end
  m = m(:);
  if length(m) > 1
    l = l(:).';
  end
  idx = m .* (2*lmax+1 - m) / 2 + l + 1;
  idx = idx(:);
end
