classdef healmex < matlab.mixin.CustomDisplay
  methods (Static)
    ipix = nest2ring(nside, ipix)
    ipix = ring2nest(nside, ipix)
    npix = nside2npix(nside)
    nside = npix2nside(npix)

    vec = pix2vec(nside, order, ipix)
    [z, phi] = pix2zphi(nside, order, ipix)
    [theta, phi] = pix2ang(nside, order, ipix)
    ipix = vec2pix(nside, order, vec)
    ipix = zphi2pix(nside, order, z, phi)
    ipix = ang2pix(nside, order, theta, phi)

    alms = map2alm(map, order, lmax, mmax, nside, niter)
    alms = map2alm_pure(map, wmap, order, lmax, mmax, nside, niter, pureE)
    map = alm2map(alms, nside, order, lmax, mmax)
    alms = hpx_map2alm(nside, order, map, lmax, mmax, rwghts, iter)
    map = hpx_alm2map(lmax, mmax, alms, nside, order)
    [almsT,almsG,almsC] = hpx_map2alm_pol(nside, order, mapT, mapQ, mapU, lmax, mmax, rwghts, iter)
    [almsG,almsC] = hpx_map2alm_pure(nside, order, mapQ, mapU, mapW, lmax, mmax, rwghts, iter, pureE)
    [mapT,mapQ,mapU] = hpx_alm2map_pol(lmax, mmax, almsT, almsG, almsC, nside, order)
    [mapQ,mapU] = hpx_alm2map_polonly(lmax, mmax, almsG, almsC, nside, order)
	
    [mapQ,mapU] = hpx_smoothing(mapQ, mapU, fle, flb, order, lmax, mmax, mmin, nside, rwghts, niter)
    map = smoothing(map, fl, mask, order, lmax, mmax, mmin, nside, niter)

    [lmax, mmax] = alm_getlmmax(alms, lmax, mmax)
    [l,m] = alm_getlm(lmax,idx)
    nel = alm_getn(lmax, mmax)
    idx = alm_getidx(lmax, l, m)

    cl = alm2cl(alms1, alms2, lmax, mmax)
    alms = almxfl(lmax, mmax, alms, fl)
    alms = rotate_alm(transform, alms, lmax, mmax)
    alms = hpx_rotate_alm(transform, alms, lmax, mmax)
    [almsT,almsG,almsC] = hpx_rotate_alm_pol(transform, almsT, almsG, almsC, lmax, mmax)

    amap = apodize_mask(map, radius, order)
    amap = hpx_apodize_mask(nside, order, map, radius)
    amap = shrink_mask(map, radius, order)
    amap = hpx_shrink_mask(nside, order, map, radius)
    amap = smooth_mask(map, radius, order)
    amap = hpx_smooth_mask(nside, order, map, radius, rwghts)
	
	rwghts = scan_rings_observed(map)
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
