function [z, phi] = pix2zphi(nside, ipix, opt)
% [z, phi] = pix2zphi(nside, ipix, ...)
%
% INPUTS
%   nside       The HEALPix Nside parameter.
%   ipix        Pixel indices.
%
% KEY-VALUE PAIRS
%   'nest'      Defaults to false. If true, `ipix` are NESTED ordering pixels,
%               otherwise assumes RING ordering.
%
% OUTPUT
%   z           The cosine of the colatitude (-1 <= z <= 1).
%   phi         The azimuth phi in radians (0 <= phi < 2*pi).
%
% EXAMPLE
%   [z, phi] = healmex.pix2zphi(512, (1000:2000)');

  arguments
    nside     (1,1) {mustBeNumeric}
    ipix            {mustBeNumeric}
    opt.nest  (1,1) logical = false
  end
  [z, phi] = libhealmex(int64(12), ...
      int64(nside), logical(opt.nest), int64(ipix));
end

