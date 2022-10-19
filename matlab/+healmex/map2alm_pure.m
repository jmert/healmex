function alms = map2alm_pure(map, wmap, order, lmax, mmax, nside, niter, pureE)
% alms = map2alm_pure(map, wmap, order, lmax, mmax, nside, niter)
%
% INPUTS
%   map     An Nx1 or Nx3 matrix of map pixel values. Column 1 is assumed to
%           be the intensity map, and columns 2 and 3 are assumed to be the
%           Stokes Q and U fields, respectively.
%
%   wmap    An Nx1 or Nx2 matrix of map pixel values. Column 1 is assumed to
%           be the intensity map weight, and column 2 is assumed to be the weight 
%           of the Stokes Q and U fields.
%
%   lmax    The maximum degree harmonic coefficient to compute. Defaults to
%           3*nside-1 if not provided.
%
%   mmax    The maximum order harmonic coefficient to compute. Defaults to
%           lmax if not provided.
%
%   nside   HEALPix nside of the map. Inferred from length of dimension one
%           of map if not provided.
%
%   niter   Number of iterations to perform in convergence to spherical
%           harmonic coefficients. Defaults to 1.
%
% OUTPUTS
%   alms    The spherical harmonic coefficients of map. If map is Nx1, then
%           alms is an Mx1 vector of scalar (non-spin) coefficients (i.e.
%           temperature). If map is Nx3, the column 1 is the scalar (non-spin)
%           coefficients, and columsn 2 and 3 are the gradient and curl,
%           respectively, spin-2 transforms.
%
% EXAMPLE
%

  if ~exist('order', 'var') || isempty(order)
    order = 'RING';
  end
  if ~exist('niter', 'var') || isempty(niter)
    niter = 1;
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

  if ~exist('pureE', 'var') || isempty(pureE)
    pureE = false;
  end


  % TODO: Allow real ring weights.
  rwghts = ones(4 * nside - 1, 1);

  if size(map, 2) == 1
    alms = healmex.hpx_map2alm(nside, order, map.*wmap, ...
        lmax, mmax, rwghts, niter);
  elseif size(map, 2) == 3
    alms(:,1) = healmex.hpx_map2alm(nside, order, map(:,1).*wmap(:,1), ...
        lmax, mmax, rwghts, niter);
    [alms(:,2),alms(:,3)] = healmex.hpx_map2alm_pure(...
        nside, order, map(:,2), map(:,3), wmap(:,2), ...
        lmax, mmax, rwghts, niter, pureE);
  else
    error('map: Expected size 1 or 3 in second dimension, got %d', ...
        size(map, 2));
  end
end
