function ipix = nest2ring(nside, ipix)
% ipix = nest2ring(nside, ipix)
%
% Converts the NESTED-ordered pixel ipix to its RING-order counterpart for
% an Nside = nside HEALPix map.
%

  ipix = libhealmex(int64(1), ...
      int64(nside), int64(ipix));
end

