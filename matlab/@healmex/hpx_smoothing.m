function [mapQ,mapU] = hpx_smoothing(mapQ, mapU, fle, flb, order, lmax, mmax, mmin, nside, rwghts, niter)
  % amap = smooth_mask(nside, order, map, radius)
  
  [mapQ,mapU] = libhealmex(int64(72), int64(nside), char(order), double(mapQ), double(mapU), double(fle), double(flb), int32(lmax), int32(mmax), int32(mmin), double(rwghts), int32(niter));
end
