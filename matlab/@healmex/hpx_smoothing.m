function map = hpx_smoothing(map, fl, order, lmax, mmax, mmin, nside, rwghts, niter)
  % amap = smooth_mask(nside, order, map, radius)
  
  map = libhealmex(int64(73), int64(nside), char(order), double(map), double(fl), int32(lmax), int32(mmax), int32(mmin), double(rwghts), int32(niter));
end
