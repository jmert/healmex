function ipix = nest2ring(nside, ipix)
% ipix = nest2ring(nside, ipix)
%
% Converts the NESTED-ordered pixel ipix to its RING-order counterpart for
% an Nside = nside HEALPix map.
%

  ipix = libhealmex(healmex.id_nest2ring, ...
      int64(nside), int64(ipix));
end

