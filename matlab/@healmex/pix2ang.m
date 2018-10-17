function [theta,phi] = pix2ang(nside, ipix)
% [theta,phi] = pix2ang(nside, ipix)
%
% Calculates HEALPix pixel center locations for pixel indices ipix in an
% Nside = nside map, returning the colatitude theta and azimuth phi spherical
% locations.
%

  [theta, phi] = libhealmex(healmex.id_pix2ang, ...
      int64(nside), int64(ipix));
end

