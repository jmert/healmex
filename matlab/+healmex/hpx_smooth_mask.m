function amap = hpx_smooth_mask(nside, order, map, radius, rwghts)
  amap = libhealmex(int64(70), int64(nside), char(order), double(map), double(radius), double(rwghts));
end
