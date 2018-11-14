function [theta,phi] = pix2ang(nside, order, ipix)
% [theta,phi] = pix2ang(nside, order, ipix)
%
% Calculates HEALPix pixel center locations for pixel indices ipix in an
% Nside = nside map with ordering scheme order, returning the colatitude theta
% and azimuth phi spherical coordinates in radians. order may be 'RING' or
% 'NESTED'.

  if ~exist('order','var') || isempty(order)
    order = 'RING';
  end

  [theta, phi] = libhealmex(int64(13), ...
      int64(nside), char(order), int64(ipix));
end

