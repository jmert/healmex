function alms = map2alm_iter(nside, order, map, lmax, mmax, rwghts, iter)
% alms = map2alm_iter(nside, order, map, lmax, mmax, rwghts, iter)
%
% Computes the spherical harmonic transform of map and returns the harmonic
% coefficients alms.

  alms = libhealmex(int64(53), ...
      int64(nside), char(order), double(map), int32(lmax), int32(mmax), ...
      double(rwghts), int32(iter));
end
