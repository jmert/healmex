function [x, y, f] = pix2xyf(nside, ipix, varargin)
% [x, y, f] = pix2xyf(nside, ipix, varargin)
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
%   x
%   y
%   f
%
% EXAMPLE
%   [x, y, f] = healmex.pix2xyf(512, (1000:2000)');

  p = inputParser();
  addParameter(p, 'nest', false, @islogical);
  parse(p, varargin{:});
  opt = p.Results;

  if opt.nest
    order = 'NESTED';
  else
    order = 'RING';
  end
  [x, y, f] = libhealmex(int64(17), ...
      int64(nside), char(order), int64(ipix));
end
