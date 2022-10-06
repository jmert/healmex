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
    nel = size(alms, 1);
    if isempty(mmax)
      % Translated from healpy's _sphtools.pyx source
      %function lmax=get_lmax1(nel)
      lmax = (sqrt(1 + 8*nel) - 3) / 2;
      if lmax ~= floor(lmax)
        throwAsCaller(MException('healmex:alm_getlmmax:sizeError', 'Could not infer lmax'));
      end
      %end
      mmax = lmax;
    else
      % Translated from healpy's _sphtools.pyx source
      %function lmax=get_lmax2(nel, mmax)
      lmax = (2*nel + mmax^2 - mmax - 2) / (2*mmax + 2);
      if lmax ~= floor(lmax)
        throwAsCaller(MException('healmex:alm_getlmmax:sizeError', 'Could not infer lmax'));
      end
      %end
    end
  elseif isempty(mmax)
    mmax = lmax;
  end
end
