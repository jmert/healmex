function ring = pix2vec(nside, ipix, opt)
% ring = pix2vec(nside, ipix, ...)
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
%   ring        Containing ring (1 <= ring <= 4*nside-1).
%
% EXAMPLE
%   ring = healmex.pix2ring(512, 0:4*512-1);

  arguments
    nside     (1,1) {mustBeNumeric}
    ipix            {mustBeNumeric}
    opt.nest  (1,1) logical = false
  end
  ring = libhealmex(int64(10), ...
      int64(nside), logical(opt.nest), int64(ipix));
end
