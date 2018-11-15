function alms = almxfl(lmax, mmax, alms, fl)
% alms = almxfl(lmax, mmax, alms, fl)
%
% Multiplies a set of alms by a vector fl.
%
% INPUTS
%   lmax    The maximum degree of the harmonic coefficients, alms.
%
%   mmax    The maximum order of the harmonic coefficients, alms.
%
%   alms    Vector of harmonic coefficients.
%
%   fl      Vector of factors f_l by which to multiply a_lm.
%
% OUTPUTS
%   alms    Modified alms.

  alms = libhealmex(int64(62), ...
      int32(lmax), int32(mmax), complex(double(alms)), complex(double(fl)));
end
