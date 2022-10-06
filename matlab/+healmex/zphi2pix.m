function ipix = zphi2pix(nside, z, phi, opt)
% ipix = zphi2pix(nside, z, phi, ...)
%
% INPUTS
%   nside       The HEALPix Nside parameter.
%   z           The cosine of the colatitude (-1 <= z <= 1).
%   phi         The azimuth phi in radians (0 <= phi < 2*pi).
%
% KEY-VALUE PAIRS
%   'nest'      Defaults to false. If true, `ipix` are NESTED ordering pixels,
%               otherwise assumes RING ordering.
%
% OUTPUT
%   ipix        Pixel indices.
%
% EXMAPLE

  arguments
    nside     (1,1) {mustBeNumeric}
    z               {mustBeNumeric}
    phi             {mustBeNumeric}
    opt.nest  (1,1) logical = false
  end
  ipix = libhealmex(int64(15), ...
      int64(nside), logical(opt.nest), double(z), double(phi));
end

