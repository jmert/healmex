function [almsG,almsC] = map2alm_pure(nside, order, mapQ, mapU, mapW, lmax, mmax, rwghts, iter, pureE)
% [almsT,almsG,almsC] = map2alm_pure(nside, order, mapQ, mapU, mapW, lmax, mmax, rwghts, iter)
%
% Computes the spherical harmonic transform of map and returns the harmonic
% coefficients alms.

  [almsG,almsC] = libhealmex(int64(67), ...
      int64(nside), char(order), double(mapQ), double(mapU), double(mapW),  ...
      int32(lmax), int32(mmax), double(rwghts), int32(iter), boolean(pureE));
end