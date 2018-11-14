function alms = map2alm_iter(nside, order, map, lmax, mmax, rwghts, iter)
% alms = map2alm_iter(nside, order, map, lmax, mmax, rwghts, iter)
%
% Computes the spherical harmonic transform of map and returns the harmonic
% coefficients alms.

  if ~exist('iter','var') || isempty(iter)
    iter = int32(0);
  end
  if ~exist('rwghts','var') || isempty(rwghts)
    rwghts = ones(4*nside - 1, 1);
  end
  if ~exist('lmax', 'var') || isempty(lmax)
    lmax = 3 * nside - 1;
  end
  if ~exist('mmax', 'var') || isempty(mmax)
    mmax = lmax;
  end

  alms = libhealmex(int64(53), ...
      int64(nside), char(order), double(map), int32(lmax), int32(mmax), ...
      double(rwghts), int32(iter));
end