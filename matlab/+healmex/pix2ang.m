function [theta, phi] = pix2ang(nside, ipix, opt)
% [theta, phi] = pix2ang(nside, ipix, ...)
%
% INPUTS
%   nside       The HEALPix Nside parameter.
%   ipix        Pixel indices.
%
% KEY-VALUE PAIRS
%   'latlon'    Defaults to false. If true, instead returns the latitude and
%               longitude coordinates [lat,lon] in degrees.
%   'nest'      Defaults to false. If true, `ipix` are NESTED ordering pixels,
%               otherwise assumes RING ordering.
%
% OUTPUT
%   theta       If latlon==true, the latitude in degrees (-90 <= lat <= 90),
%               otherwise the colatitude in radians (0 <= theta <= pi).
%   phi         If latlon==true, the longitude in degrees (0 <= lon < 360),
%               otherwise the azimuth in radians (0 <= phi < 2*pi).
%
% EXAMPLE
%   [lat,lon] = healmex.pix2ang(512, (1000:2000)', 'latlon', true);

  arguments
    nside       (1,1) {mustBeNumeric}
    ipix              {mustBeNumeric}
    opt.latlon  (1,1) logical = false
    opt.nest    (1,1) logical = false
  end

  [theta, phi] = libhealmex(int64(13), ...
      int64(nside), logical(opt.nest), int64(ipix));
  if opt.latlon
    theta = 90 - rad2deg(theta);
    phi = rad2deg(phi);
  end
end

