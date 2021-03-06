classdef healmex < matlab.mixin.CustomDisplay
  methods (Static)
    ipix = nest2ring(nside, ipix)
    ipix = ring2nest(nside, ipix)
    npix = nside2npix(nside)
    nside = npix2nside(npix)

    vec = pix2vec(nside, order, ipix)
    [z, phi] = pix2zphi(nside, order, ipix)
    [theta, phi] = pix2ang(nside, order, ipix)
    [x, y, f] = pix2xyf(nside, order, ipix)
    ipix = vec2pix(nside, order, vec)
    ipix = zphi2pix(nside, order, z, phi)
    ipix = ang2pix(nside, order, theta, phi)
    ipix = xyf2pix(nside, order, x, y, f)

    alms = map2alm(map, order, lmax, mmax, nside, niter)
    map = alm2map(alms, nside, order, lmax, mmax)
    alms = hpx_map2alm(nside, order, map, lmax, mmax, rwghts, iter)
    [almsT,almsG,almsC] = hpx_map2alm_pol(nside, order, mapT, mapQ, mapU, lmax, mmax, rwghts, iter)

    [lmax, mmax] = alm_getlmmax(alms, lmax, mmax)
    nel = alm_getn(lmax, mmax)
    idx = alm_getidx(lmax, l, m)

    cl = alm2cl(alms1, alms2, lmax, mmax)
    alms = almxfl(lmax, mmax, alms, fl)
    alms = rotate_alm(transform, alms, lmax, mmax)
  end

  methods (Access = protected)
    function str=getHeader(self)
      str = [class(self) ' with methods:'];
    end
    function P=getPropertyGroups(self)
      persistent pg
      if isempty(pg)
        mt = setdiff(methods(self), class(self));
        pg = struct();
        for ii = 1:length(mt)
          docs = help([class(self) '.' mt{ii}]);
          % char(10) == '\n' without needing to be interpreted
          pg.(mt{ii}) = strtrim(strtok(docs, char(10)));
        end
      end
      P = matlab.mixin.util.PropertyGroup(pg);
    end
  end
end
