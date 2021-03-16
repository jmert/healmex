classdef healmex < matlab.mixin.CustomDisplay
  properties (Constant)
    UNSEEN = single(-1.6375e30)
  end

  methods (Static)
    ipix  = nest2ring(nside, ipix)
    ipix  = ring2nest(nside, ipix)
    npix  = nside2npix(nside)
    nside = npix2nside(npix)

    [theta, phi] = pix2ang(nside, ipix, varargin)
    [x, y, z]    = pix2vec(nside, ipix, varargin)
    [x, y, f]    = pix2xyf(nside, ipix, varargin)
    [z, phi]     = pix2zphi(nside, ipix, varargin)

    ipix = ang2pix(nside, theta, phi, varargin)
    ipix = vec2pix(nside, x, y, z, varargin)
    ipix = xyf2pix(nside, x, y, f, varargin)
    ipix = zphi2pix(nside, z, phi, varargin)

    alms = map2alm(map, varargin)
    maps = alm2map(alms, nside, varargin)
    maps = alm2map_der1(alms, nside, varargin)

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
