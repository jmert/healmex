function map = hpx_alm2map(lmax, mmax, alms, nside, order)
% map = hpx_alm2map(lmax, mmax, alms, nside, order)
%
% Synthesizes a map at Nside = nside and pixel ordering order from the set of
% spherical harmonic coefficients alms.

  map = libhealmex(int64(55), ...
      int32(lmax), int32(mmax), complex(double(alms)), ...
      int64(nside), char(order));
end
