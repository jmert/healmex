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

  nstokes=size(map, 2);
  
  if nstokes>= 2
    if ~exist('mask', 'var') || isempty(mask)
      [map(:,nstokes-1),map(:,nstokes)] = healmex.hpx_smoothing_pol(map(:,nstokes-1), -1.*map(:,nstokes), fle, flb, order, lmax, mmax, mmin, nside, rwghts, niter);
	else
      [map(:,nstokes-1),map(:,nstokes)] = healmex.hpx_smoothing_pol(map(:,nstokes-1).*mask, -1.*map(:,nstokes).*mask, fle, flb, order, lmax, mmax, mmin, nside, rwghts, niter);
	end
    map(:,nstokes)=-1.*map(:,nstokes);
  else if nstokes == 3
    if ~exist('mask', 'var') || isempty(mask)
      map(:,1) = healmex.hpx_smoothing(map(:,1), fl, order, lmax, mmax, mmin, nside, rwghts, niter);
	else
      map(:,1) = healmex.hpx_smoothing(map(:,1).*mask, fl, order, lmax, mmax, mmin, nside, rwghts, niter);
	end
  else
    error('map: Expected size 2 or 3 in second dimension, got %d', nstokes);
  end
end
