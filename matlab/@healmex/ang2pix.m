function ipix = ang2pix(nside, order, theta, phi)
% ipix = ang2pix(nside, order, theta, phi)
%
% Calculates HEALPix pixel indices ipix (for an Nside = nside map with ordering
% scheme order) which contain the spherical coordinate points (theta,phi),
% where theta and phi are colatitude/azimuth angle in radians. order may be
% 'RING' or 'NESTED'.

  if ~exist('order','var') || isempty(order)
    order = 'RING';
  end

  ipix = libhealmex(int64(16), ...
      int64(nside), char(order), double(theta), double(phi));
end

