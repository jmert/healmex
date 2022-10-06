function ipix = ang2pix(nside, theta, phi, opt)
% ipix = ang2pix(nside, theta, phi, ...)
%
% INPUTS
%   nside       The HEALPix Nside parameter.
%   theta       If latlon==true, the latitude in degrees (-90 <= lat <= 90),
%               otherwise the colatitude in radians (0 <= theta <= pi).
%   phi         If latlon==true, the longitude in degrees (0 <= lon < 360),
%               otherwise the azimuth in radians (0 <= phi < 2*pi).
%
% KEY-VALUE PAIRS
%   'latlon'    Defaults to false. If true, instead interprets theta and phi
%               as latitude and longitude coordinates in degrees.
%   'nest'      Defaults to false. If true, `ipix` are NESTED ordering pixels,
%               otherwise assumes RING ordering.
%
% OUTPUT
%   ipix        Pixel indices.
%

  arguments
    nside       (1,1) {mustBeNumeric}
    theta             {mustBeNumeric}
    phi               {mustBeNumeric}
    opt.latlon  (1,1) logical = false
    opt.nest    (1,1) logical = false
  end
  if opt.latlon
    theta = deg2rad(90 - theta);
    phi = deg2rad(phi);
  end
  ipix = libhealmex(int64(16), ...
      int64(nside), logical(opt.nest), double(theta), double(phi));
end

