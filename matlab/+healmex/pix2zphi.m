function [z, phi] = pix2zphi(nside, ipix, varargin)
% [z, phi] = pix2zphi(nside, ipix, varargin)
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

  p = inputParser();
  addParameter(p, 'nest', false, @islogical);
  parse(p, varargin{:});
  opt = p.Results;

  if opt.nest
    order = 'NESTED';
  else
    order = 'RING';
  end
  [z, phi] = libhealmex(int64(12), ...
      int64(nside), char(order), int64(ipix));
end

