function ipix = vec2pix(nside, x, y, z, varargin)
% ipix = vec2pix(nside, order, x, y, z, varargin)
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

  p = inputParser();
  addParameter(p, 'nest', false, @islogical);
  parse(p, varargin{:});
  opt = p.Results;

  if opt.nest
    order = 'NESTED';
  else
    order = 'RING';
  end
  ipix = libhealmex(int64(14), ...
      int64(nside), char(order), double(x), double(y), double(z));
end

