function idx = alm_getidx(lmax, l, m)
% idx = alm_getidx(lmax, l, m)
%
% Returns the index of the alm coefficient (l,m) within a vector which includes
% coefficients up to lmax.

<<<<<<< HEAD:matlab/@healmex/alm_getidx.m
  if ~isscalar(lmax)
    error('lmax must be a scalar');
  end
  
=======
  arguments
    lmax  (1,1) {mustBeInteger}
    l     (:,:) {mustBeNumeric}
    m     (:,:) {mustBeNumeric}
  end

>>>>>>> master:matlab/+healmex/alm_getidx.m
  idx = m .* (2*lmax+1 - m) / 2 + l + 1;
end
