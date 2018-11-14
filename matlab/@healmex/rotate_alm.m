function alms=rotate_alm(transform, alms, lmax, mmax)
% alms=rotate_alm(transform, alms, lmax, mmax)
%
% Performs coordinate transformation (e.g. rotations on the sphere) in
% harmonic space.
%
% INPUTS
%   transform   An integer giving one of the following coordinate transforms.
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
%   alms        A vector of spherical harmonic coefficients.
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
  [lmax, mmax] = healmex.alm_getlmmax(alms, lmax, mmax);

  if transform < 1 || transform > 12
    error('Unknown coordinate transform: %d', transform)
  end

  alms = libhealmex(int64(65), ...
      int32(transform), int32(lmax), int32(mmax), complex(double(alms)));
end
