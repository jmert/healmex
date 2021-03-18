function [theta, phi] = pix2ang(nside, ipix, opt)
% [theta, phi] = pix2ang(nside, ipix, ...)
%
% INPUTS
%   nside       The HEALPix Nside parameter.
%   ipix        Pixel indices.
%
% KEY-VALUE PAIRS
%   'lonlat'    Defaults to false. If true, instead returns the longitude and
%               latitude coordinates [lon,lat] in degrees.
%   'nest'      Defaults to false. If true, `ipix` are NESTED ordering pixels,
%               otherwise assumes RING ordering.
%
% OUTPUT
%   theta       If lonlat==true, the latitude in degrees (-90 <= lat <= 90),
%               otherwise the colatitude in radians (0 <= theta <= pi).
%   phi         If lonlat==true, the longitude in degrees (0 <= lon < 360),
%               otherwise the azimuth in radians (0 <= phi < 2*pi).
%
% EXAMPLE
%   [lon,lat] = healmex.pix2ang(512, (1000:2000)', 'lonlat', true);

  arguments
    nside       (1,1) {mustBeNumeric}
    ipix              {mustBeNumeric}
    opt.lonlat  (1,1) logical = false
    opt.nest    (1,1) logical = false
  end

  if opt.nest
    order = 'NESTED';
  else
    order = 'RING';
  end

  [theta, phi] = libhealmex(int64(13), ...
      int64(nside), char(order), int64(ipix));
  if opt.lonlat
    theta = 90 - rad2deg(theta);
    phi = rad2deg(phi);
  end
end

