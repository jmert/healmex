function amap = smooth_mask(nside, order, map, radius)
% amap = smooth_mask(nside, order, map, radius)

  amap = libhealmex(int64(69), ...
      int64(nside), char(order), double(map), double(radius));
end
