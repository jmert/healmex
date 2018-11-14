function [almsT,almsG,almsC]=rotate_alm_pol(transform, almsT, almsG, almsC, lmax, mmax)
% [almsT,almsG,almsC]=rotate_alm_pol(transform, almsT, almsG, almsC, lmax, mmax)
%
% Performs coordinate transformation (e.g. rotations on the sphere) in
% harmonic space. Exactly same as rotate_alm(), but performs rotations on
% three sets of alms simultaneously.
%
% INPUTS
%   transform
%
%       1 = Equatorial (2000) -> Galactic   (2000)
%       2 = Galactic   (2000) -> Equatorial (2000)
%       3 = Equatorial (2000) -> Ecliptic   (2000)
%       4 = Ecliptic   (2000) -> Equatorial (2000)
%       5 = Ecliptic   (2000) -> Galactic   (2000)
%       6 = Galactic   (2000) -> Ecliptic   (2000)
%       7 = Equatorial (1950) -> Galactic   (1950)
%       8 = Galactic   (1950) -> Equatorial (1950)
%       9 = Equatorial (1950) -> Ecliptic   (1950)
%      10 = Ecliptic   (1950) -> Equatorial (1950)
%      11 = Ecliptic   (1950) -> Galactic   (1950)
%      12 = Galactic   (1950) -> Ecliptic   (1950)
%
%   almsT       A vector of spherical harmonic coefficients.
%   almsG       A vector of spherical harmonic coefficients.
%   almsC       A vector of spherical harmonic coefficients.
%   lmax        Maximum degree of spherical harmonics. Optional if inferrable
%               by alm_getlmmax().
%   mmax        Maximum order of spherical harmonics. Optional if inferrable
%               by alm_getlmmax().

  if ~exist('lmax', 'var')
    lmax = [];
  end
  if ~exist('mmax', 'var')
    mmax = [];
  end
  [lmax, mmax] = healmex.alm_getlmmax(almsT, lmax, mmax);
  if numel(almsT) ~= numel(almsG)
    error('Mismatched sizes in almsT and almsG');
  end
  if numel(almsT) ~= numel(almsC)
    error('Mismatched sizes in almsT and almsC');
  end

  if iscell(transform) && length(transform) == 2
  end

  [almsT,almsG,almsC] = libhealmex(healmex.id_rotate_alm_pol, ...
      int32(transform), int32(lmax), int32(mmax), ...
      complex(double(almsT)), complex(double(almsG)), complex(double(almsC)));
end
