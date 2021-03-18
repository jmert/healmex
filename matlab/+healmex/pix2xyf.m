function [x, y, f] = pix2xyf(nside, ipix, opt)
% [x, y, f] = pix2xyf(nside, ipix, ...)
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
  [x, y, f] = libhealmex(int64(17), ...
      int64(nside), char(order), int64(ipix));
end
