function [mapT,mapQ,mapU] = alm2map_pol(lmax, mmax, almsT, almsG, almsC, nside, order)
% [mapT,mapQ,mapU] = alm2map_pol(lmax, mmax, almsT, almsG, almsC, nside, order)
%
% Synthesizes a map at Nside = nside and pixel ordering order from the set of
% spherical harmonic coefficients alms.

  if ~exist('order', 'var') || isempty(order)
    order = 'RING';
  end

  [mapT, mapQ, mapU] = libhealmex(healmex.id_alm2map_pol, ...
      int32(lmax), int32(mmax), double(almsT), double(almsG), double(almsC), ...
      int64(nside), char(order));
end
