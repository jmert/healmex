function ipix = ring2nest(nside, ipix)
% ipix = ring2nest(nside, ipix)
%
% INPUTS
%   nside     HEALPix Nside resolution factor.
%   ipix      Pixel indices in RING ordering.
%
% OUTPUTS
%   ipix      Pixel indices in NESTED ordering.
%
% EXAMPLE
%   ringpix = healmex.ring2nest(512, nestpix);

  ipix = libhealmex(int64(2), ...
      int64(nside), int64(ipix));
end
