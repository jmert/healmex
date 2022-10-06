function amap = apodize_mask(map, radius, order)
% amap = apodize_mask(map, radius, order)
%
% INPUTS
%   map     An Nx1 matrix of mask pixel values.
%
%   radius  Radius of edge taper, using the analytical function 0.5 + (0.5* cos(delta))
%           (see Grain et al. 2009)
%
% OUTPUTS
%   amap    Apodized mask
%
% EXAMPLE
%

  if ~exist('order', 'var') || isempty(order)
    order = 'RING';
  end
  
  nside = healmex.npix2nside(size(map, 1));
  
  amap=healmex.hpx_apodize_mask(nside,order,map,radius);
end
