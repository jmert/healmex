function r = pix2vec(nside, order, ipix)
% r = pix2vec(nside, order, ipix)
%
% Calculates HEALPix pixel center locations for pixel indices ipix in an
% Nside = nside map with ordering scheme order, returning the unit vector r
% pointing to the pixel center as an N-by-3 matrix of x, y, and z
% Cartesian coordinates. order may be 'RING' or 'NESTED'.

  r = libhealmex(int64(11), ...
      int64(nside), char(order), int64(ipix));
end

