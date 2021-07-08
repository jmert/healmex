function ipix = ring2nest(nside, ipix)
% ipix = ring2nest(nside, ipix)
%
% Converts the RING-ordered pixel ipix to its NESTED-order counterpart for
% an Nside = nside HEALPix map.
%

  ipix = libhealmex(int64(2), ...
      int64(nside), int64(ipix));
end

