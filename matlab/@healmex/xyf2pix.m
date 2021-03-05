function ipix = xyf2pix(nside, order, x, y, f)
% ipix = pix2xyf(nside, order, x, y, f)
%
% Calculates HEALPix pixel indices ipix in an Nside = nside map with ordering
% scheme order corresponding to locations identified by the triples (x,y,face).
% order may be 'RING' or 'NESTED'.

  if ~exist('order','var') || isempty(order)
    order = 'RING';
  end

  ipix = libhealmex(int64(18), ...
      int64(nside), char(order), int32(x), int32(y), int32(f));
end
