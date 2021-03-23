function nbrpix = neighbors(nside, ipix, opt)
% nbrpix = neightbor(nside, ipix, ...)
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
%   nbrpix      Neighboring SW, W, NW, N, NE, E, SE, and S pixels as an
%               8-by-length(ipix) matrix. If a neighbor does not exist, the
%               corresponding pixel value will be -1.
%
% EXAMPLE
%   ring = healmex.neighbors(1, 4);

  arguments
    nside     (1,1) {mustBeNumeric}
    ipix            {mustBeNumeric}
    opt.nest  (1,1) logical = false
  end
  nbrpix = libhealmex(int64(20), ...
      int64(nside), logical(opt.nest), int64(ipix));
end
