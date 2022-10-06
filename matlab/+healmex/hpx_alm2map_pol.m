function [mapT,mapQ,mapU] = hpx_alm2map_pol(lmax, mmax, almsT, almsG, almsC, nside, order)
% [mapT,mapQ,mapU] = hpx_alm2map_pol(lmax, mmax, almsT, almsG, almsC, nside, order)
%
% Synthesizes a map at Nside = nside and pixel ordering order from the set of
% spherical harmonic coefficients alms.

  [mapT, mapQ, mapU] = libhealmex(int64(56), ...
      int32(lmax), int32(mmax), ...
      complex(double(almsT)), complex(double(almsG)), complex(double(almsC)), ...
      int64(nside), char(order));
end
