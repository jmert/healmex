function amap = apodize_mask(map, radius, nest)
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

  if ~exist('nest', 'var') || isempty(nest)
    nest = false;
  end
  
  nside = healmex.npix2nside(size(map, 1));
  
  amap = libhealmex(int64(68), int64(nside), logical(nest), double(map), double(radius));
end