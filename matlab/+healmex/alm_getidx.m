function idx = alm_getidx(lmax, l, m)
% idx = alm_getidx(lmax, l, m)
%
% Returns the index of the alm coefficient (l,m) within a vector which includes
% coefficients up to lmax.

  arguments
    lmax  (1,1) {mustBeInteger}
    l     (:,:) {mustBeNumeric}
    m     (:,:) {mustBeNumeric}
  end

  idx = m .* (2*lmax+1 - m) / 2 + l + 1;
end
