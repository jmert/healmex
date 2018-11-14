function [z,phi] = pix2zphi(nside, order, ipix)
% [z,phi] = pix2zphi(nside, order, ipix)
%
% Calculates HEALPix pixel center locations for pixel indices ipix in an
% Nside = nside map with ordering scheme order, returning the cosine of the
% colatitude theta (i.e the Cartesian z coordinate of a unit vector) and
% azimuth angle phi (in radians). order may be 'RING' or 'NESTED'.

  if ~exist('order','var') || isempty(order)
    order = 'RING';
  end

  [z, phi] = libhealmex(int64(12), ...
      int64(nside), char(order), int64(ipix));
end

