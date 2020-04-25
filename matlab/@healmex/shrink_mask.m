function amap = shrink_mask(map, radius, order)
% amap = shrink_mask(map, radius, order)
%
% INPUTS
%   map     An Nx1 matrix of mask pixel values.
%
%   radius  Size of perimeter in degree
%           
%
% OUTPUTS
%   amap    Mask with perimeter cut out
%
% EXAMPLE
%

  if ~exist('order', 'var') || isempty(order)
    order = 'RING';
  end
  
  nside = healmex.npix2nside(size(map, 1));
  
  amap=healmex.hpx_shrink_mask(nside,order,map,radius);
end
