function [lmax, mmax] = alm_getlmmax(alms, lmax, mmax)
% [lmax, mmax] = alm_getlmmax(alms, lmax, mmax)
%
% Infers the lmax and/or mmax from the alms vector (or length of the
% first dimension if alms is a matrix).

  arguments
    alms (:,:)
    lmax       {mustBeNumeric,mustBeScalarOrEmpty} = []
    mmax       {mustBeNumeric,mustBeScalarOrEmpty} = []
  end

  if isempty(lmax)
    if isempty(mmax)
      lmax = get_lmax1(size(alms,1));
      mmax = lmax;
    else
      lmax = get_lmax2(size(alms,1), mmax);
    end
  elseif isempty(mmax)
    mmax = lmax;
  end

end
%
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
