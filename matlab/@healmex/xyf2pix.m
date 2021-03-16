function ipix = xyf2pix(nside, x, y, f, varargin)
% ipix = xyf2pix(nside, x, y, f, varargin)
%
% INPUTS
%   nside       The HEALPix Nside parameter.
%   x
%   y
%   f
%
% KEY-VALUE PAIRS
%   'nest'      Defaults to false. If true, `ipix` are NESTED ordering pixels,
%               otherwise assumes RING ordering.
%
% OUTPUT
%   ipix        Pixel indices.
%
% EXAMPLE
%   ipix = healmex.pix2xyf(512, x, y, f);

  p = inputParser();
  addParameter(p, 'nest', false, @islogical);
  parse(p, varargin{:});
  opt = p.Results;

  if opt.nest
    order = 'NESTED';
  else
    order = 'RING';
  end
  ipix = libhealmex(int64(18), ...
      int64(nside), char(order), int32(x), int32(y), int32(f));
end
