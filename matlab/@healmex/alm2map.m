function map = alm2map(lmax, mmax, alms, nside, order)
% map = alm2map(lmax, mmax, alms, nside, order)
%
% Synthesizes a map at Nside = nside and pixel ordering order from the set of
% spherical harmonic coefficients alms.

  if ~exist('order', 'var') || isempty(order)
    order = 'RING';
  end

  map = libhealmex(int64(55), ...
      int32(lmax), int32(mmax), complex(double(alms)), ...
      int64(nside), char(order));
end
