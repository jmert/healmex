function [almsT,almsG,almsC] = map2alm_pol_iter(nside, order, mapT, mapQ, mapU, lmax, mmax, rwghts, iter)
% [almsT,almsG,almsC] = map2alm_pol_iter(nside, order, mapT, mapQ, mapU, lmax, mmax, rwghts, iter)
%
% Computes the spherical harmonic transform of map and returns the harmonic
% coefficients alms.

  [almsT,almsG,almsC] = libhealmex(int64(54), ...
      int64(nside), char(order), double(mapT), double(mapQ), double(mapU), ...
      int32(lmax), int32(mmax), double(rwghts), int32(iter));
end
