function cl = alm2cl(alms1, alms2, lmax, mmax)
% cl = alm2cl(alms1, alms2, lmax, mmax)
%
% Computes the angular power [cross-]spectrum from the set of harmonic
% coefficients alms1 and alms2. If alms2 is empty, then alms1 is used to
% produce the auto-spectrum. lmax and mmax are optional in cases where the
% lmax and mmax can be inferred from the length of the alms vector.

  if ~exist('lmax','var')
    lmax = [];
  end
  if ~exist('mmax','var')
    mmax = [];
  end

  [lmax,mmax] = healmex.alm_getlmmax(alms1, lmax, mmax);
  if ~exist('alms2', 'var') || isempty(alms2)
    alms2 = alms1;
  end
  if numel(alms1) ~= numel(alms2)
    error('Mismatched sizes in alms1 and alms2');
  end

  cl = libhealmex(int64(61), ...
      int32(lmax), int32(mmax), ...
      complex(double(alms1)), complex(double(alms2)));
end

