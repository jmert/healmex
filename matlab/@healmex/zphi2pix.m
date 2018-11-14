function ipix = zphi2pix(nside, order, z, phi)
% ipix = zphi2pix(nside, order, z, phi)
%
% Calculates HEALPix pixel indices ipix (for an Nside = nside map with ordering
% scheme order) which contain the coordinate points (z, phi), where z is a
% Cartesian coordinate and phi is the azimuth angle in radians. order may be
% 'RING' or 'NESTED'.

  if ~exist('order','var') || isempty(order)
    order = 'RING';
  end

  ipix = libhealmex(int64(15), ...
      int64(nside), char(order), double(z), double(phi));
end

