function [x, y, z] = pix2vec(nside, ipix, opt)
% [x, y, z] = pix2vec(nside, ipix, ...)
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
%   x           Cartesian X components of pointing unit vectors.
%   y           Cartesian Y components of pointing unit vectors.
%   z           Cartesian Z components of pointing unit vectors.
%
% EXAMPLE
%   [x, y, z] = healmex.pix2vec(512, 0:4*512-1);

  arguments
    nside     (1,1) {mustBeNumeric}
    ipix            {mustBeNumeric}
    opt.nest  (1,1) logical = false
  end

  if opt.nest
    order = 'NESTED';
  else
    order = 'RING';
  end
  [x,y,z] = libhealmex(int64(11), ...
      int64(nside), char(order), int64(ipix));
end

