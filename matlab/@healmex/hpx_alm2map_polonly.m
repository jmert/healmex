function [mapQ,mapU] = alm2map_polonly(lmax, mmax, almsG, almsC, nside, order, rwghts)
% [mapT,mapQ,mapU] = alm2map_pol(lmax, mmax, almsT, almsG, almsC, nside, order)
%
% Synthesizes a map at Nside = nside and pixel ordering order from the set of
% spherical harmonic coefficients alms.

  [mapQ, mapU] = libhealmex(int64(71), ...
      int32(lmax), int32(mmax), ...
      complex(double(almsG)), complex(double(almsC)), ...
      int64(nside), char(order), double(rwghts));
end
