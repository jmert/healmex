function nel=alm_getn(lmax, mmax)
% nel=alm_getn(lmax, mmax)
%
% Computes size of alm vector expected for maximum degree lmax, maximum order
% mmax.
%
% INPUTS
%   lmax    Must be a non-negative integer.
%   mmax    Defaults to lmax. Must be in the range 0 to lmax.

  arguments
    lmax {mustBeNonnegative}
    mmax {mustBeNonnegative} = lmax
  end
  nel = ((mmax + 1) .* (mmax + 2)) / 2 + (mmax + 1) .* (lmax - mmax);
end
