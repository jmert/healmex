function nel=alm_getn(lmax, mmax)
% nel=alm_getn(lmax, mmax)
%
% Computes size of alm vector expected for maximum degree lmax, maximum order
% mmax.
%
% INPUTS
%   lmax    Must be a non-negative integer.
%   mmax    Defaults to lmax. Must be in the range 0 to lmax.

  if lmax < 0 || fix(lmax) ~= lmax
    error('lmax must be a positive integer')
  end
  if ~exist('mmax', 'var') || isempty(mmax)
    mmax = lmax;
  end

  nel = ((mmax + 1) * (mmax + 2)) / 2 + (mmax + 1) * (lmax - mmax);
end
