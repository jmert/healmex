function ipix = nest2ring(nside, ipix)
% ipix = nest2ring(nside, ipix)
%
% INPUTS
%   nside     HEALPix Nside resolution factor.
%   ipix      Pixel indices in NESTED ordering.
%
% OUTPUTS
%   ipix      Pixel indices in RING ordering.
%
% EXAMPLE
%   nestpix = healmex.nest2ring(512, ringpix);

  ipix = libhealmex(int64(1), ...
      int64(nside), int64(ipix));
end
