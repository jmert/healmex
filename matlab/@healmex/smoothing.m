function map = smoothing(map, fl, mask, rwghts, order, lmax, mmax, mmin, nside, niter)
  if ~exist('order', 'var') || isempty(order)
    order = 'RING';
  end
  if ~exist('niter', 'var') || isempty(niter)
    niter = 3;
  end

  if ~exist('nside', 'var') || isempty(nside)
    nside = healmex.npix2nside(size(map, 1));
  end

  if ~exist('lmax', 'var') || isempty(lmax)
    lmax = 3 * nside - 1;
  end
  if ~exist('mmax', 'var') || isempty(mmax)
    mmax = lmax;
  end
  if ~exist('mmin', 'var') || isempty(mmin)
    mmin = 0;
  end
  if ~exist('rwghts', 'var') || isempty(rwghts)
    rwghts = ones(4 * nside - 1, 1);
  end

  if size(fl,2)>1
    fle=fl(:,1);
    flb=fl(:,2);
  else
    fle=fl(:,1);
    flb=fl(:,1);
  end

  if size(map, 2) == 2
    if ~exist('mask', 'var') || isempty(mask)
      [map(:,1),map(:,2)] = healmex.hpx_smoothing(map(:,1), -1.*map(:,2), fle, flb, order, lmax, mmax, mmin, nside, rwghts, niter);
	else
      [map(:,1),map(:,2)] = healmex.hpx_smoothing(map(:,1).*mask, -1.*map(:,2).*mask, fle, flb, order, lmax, mmax, mmin, nside, rwghts, niter);
	end
    map(:,2)=-1.*map(:,2);
  else
    error('map: Expected size 2 in second dimension, got %d', size(map, 2));
  end
end
