classdef healmex < matlab.mixin.CustomDisplay
  methods (Static)
    ipix = nest2ring(nside, ipix)
    ipix = ring2nest(nside, ipix)

    vec = pix2vec(nside, order, ipix)
    [z, phi] = pix2zphi(nside, order, ipix)
    [theta, phi] = pix2ang(nside, order, ipix)
    ipix = vec2pix(nside, order, vec)
    ipix = zphi2pix(nside, order, z, phi)
    ipix = ang2pix(nside, order, theta, phi)

    alms = map2alm_iter(nside, order, map, lmax, mmax, rwghts, iter)
    [almsT,almsG,almsC] = map2alm_pol_iter(nside, order, mapT, mapQ, mapU, lmax, mmax, rwghts, iter)
    map = alm2map(lmax, mmax, alms, nside, order)
    [mapT,mapQ,mapU] = alm2map_pol(lmax, mmax, almsT, almsG, almsC, nside, order)

    [lmax, mmax] = alm_getlmmax(alms, lmax, mmax);
    nel = alm_getn(lmax, mmax);

    cl = alm2cl(alms1, alms2, lmax, mmax)
    alms = rotate_alm(transform, alms, lmax, mmax);
    [almsT,almsG,almsC] = rotate_alm_pol(transform, almsT, almsG, almsC, lmax, mmax);
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
