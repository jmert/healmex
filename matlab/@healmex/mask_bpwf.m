function mcm = mask_bpwf(mask, lmax, lmax_mask, order, niter, rwghts, pureB)
% alms = map2alm(map, order, lmax, mmax, nside, niter)
%
% INPUTS
%   mask     
%           
%           
%
%   lmax    
%           
%
%   lmax_mask    
%           
%
% OUTPUTS
%   mcm    
%           
%           
%           
%           
%
% EXAMPLE
%

  if ~exist('order', 'var') || isempty(order)
    order = 'RING';
  end
  if ~exist('niter', 'var') || isempty(niter)
    niter = 1;
  end
  if ~exist('pureB', 'var') || isempty(niter)
    pureB = false;
  end

  pe1=0;
  pe2=0;
  pb1=0;
  pb2=0;

  if pureB
    pb1=1;
    pb2=1;
  end

  nside = healmex.npix2nside(size(mask, 1));

  % TODO: Allow real ring weights.
  if ~exist('rwghts', 'var') || isempty(rwghts)
    rwghts = ones(4 * nside - 1, 1);
  end

  mcm = libhealmex(int64(75), int32(lmax), int32(lmax_mask), int64(nside), char(order), double(mask), double(rwghts), int32(niter), int32(pe1), int32(pe2), int32(pb1), int32(pb2));
end