function map = alm2map(lmax, mmax, alms, nside, order)
% map = alm2map(lmax, mmax, alms, nside, order)
%
% Synthesizes a map at Nside = nside and pixel ordering order from the set of
% spherical harmonic coefficients alms.

  if ~exist('order', 'var') || isempty(order)
    order = 'RING';
  end

  map = libhealmex(healmex.id_alm2map, ...
      int32(lmax), int32(mmax), double(alms), int64(nside), char(order));
end
