function ipix = vec2pix(nside, x, y, z, opt)
% ipix = vec2pix(nside, x, y, z, ...)
%
% INPUTS
%   nside       The HEALPix Nside parameter.
%   x           Cartesian X components of pointing unit vectors.
%   y           Cartesian Y components of pointing unit vectors.
%   z           Cartesian Z components of pointing unit vectors.
%
% KEY-VALUE PAIRS
%   'nest'      Defaults to false. If true, `ipix` are NESTED ordering pixels,
%               otherwise assumes RING ordering.
%
% OUTPUT
%   ipix        Pixel indices.
%
% EXAMPLE
%   [x, y, z] = sph2cart(0.0:0.1:pi/2, 0.0:0.05:pi/4, 1);
%   ipix = healmex.pix2vec(512, x, y, z);

  arguments
    nside    (1,1) {mustBeNumeric}
    x              {mustBeNumeric}
    y              {mustBeNumeric}
    z              {mustBeNumeric}
    opt.nest (1,1) logical = false
  end
  ipix = libhealmex(int64(14), ...
      int64(nside), logical(opt.nest), double(x), double(y), double(z));
end

