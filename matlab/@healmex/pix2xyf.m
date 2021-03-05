function [x,y,f] = pix2xyf(nside, order, ipix)
% [x,y,f] = pix2xyf(nside, order, ipix)
%
% Calculates HEALPix pixel locations for pixel indices ipix in an Nside = nside
% map with ordering scheme order, returning the (x,y,face) values.
% order may be 'RING' or 'NESTED'.

  if ~exist('order','var') || isempty(order)
    order = 'RING';
  end

  [x, y, f] = libhealmex(int64(17), ...
      int64(nside), char(order), int64(ipix));
end
