function amap = apodize_mask(nside, order, map, radius)
% amap = map2alm_pol_iter(nside, order, mapT, mapQ, mapU, lmax, mmax, rwghts, iter)
%
% Computes the spherical harmonic transform of map and returns the harmonic
% coefficients alms.

  amap = libhealmex(int64(68), ...
      int64(nside), char(order), double(map), double(radius));
end
