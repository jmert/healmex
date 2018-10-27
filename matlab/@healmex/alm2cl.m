function cl = alm2cl(alms1, alms2, lmax, mmax)
% alms = map2alm_iter(alms1, alms2, lmax, mmax)
%
% Computes the spherical harmonic transform of map and returns the harmonic
% coefficients alms.

  if ~exist('lmax', 'var')
    lmax = [];
  end
  if ~exist('mmax', 'var')
    mmax = [];
  end

  if isempty(lmax)
    if isempty(mmax)
      lmax = get_lmax1(numel(alms1));
      mmax = lmax;
    else
      lmax = get_lmax2(numel(alms2), mmax);
    end
  elseif isempty(lmax)
    mmax = lmax;
  end

  if ~exist('alms2', 'var') || isempty(alms2)
    alms2 = alms1;
  end

  cl = libhealmex(healmex.id_alm2cl, ...
      int32(lmax), int32(mmax), double(alms1), double(alms2));
end

% Translated from healpy's _sphtools.pyx source
function lmax=get_lmax1(nel)
  lmax = (sqrt(1 + 8*nel) - 3) / 2;
  if lmax ~= floor(lmax)
    error('Could not infer lmax');
  end
end

% Translated from healpy's _sphtools.pyx source
function lmax=get_lmax2(nel, mmax)
  lmax = (2*nel + mmax^2 - mmax - 2) / (2*mmax + 2);
  if lmax ~= floor(lmax)
    error('Could not infer lmax');
  end
end
