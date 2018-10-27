function ipix = vec2pix(nside, order, vec)
% ipix = vec2pix(nside, order, vec)
%
% Calculates HEALPix pixel indices ipix (for an Nside = nside map with ordering
% scheme order) which contains the endpoints of the unit vectors vec, wheree
% vec is a N-by-3 matrix of Cartesian coordinates. order may be 'RING' or
% 'NESTED'.

  if ~exist('order','var') || isempty(order)
    order = 'RING';
  end

  ipix = libhealmex(healmex.id_vec2pix, ...
      int64(nside), char(order), double(vec));
end

