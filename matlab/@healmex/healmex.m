classdef healmex < matlab.mixin.CustomDisplay
  properties (Constant, Access = private)
    id_heartbeat        = int64(-1)
    id_nest2ring        = int64(1)
    id_ring2nest        = int64(2)
    id_pix2vec          = int64(11)
    id_pix2zphi         = int64(12)
    id_pix2ang          = int64(13)
    id_vec2pix          = int64(14)
    id_zphi2pix         = int64(15)
    id_ang2pix          = int64(16)
  end
  methods (Static)
    ipix = nest2ring(nside, ipix)
    ipix = ring2nest(nside, ipix)

    vec = pix2vec(nside, order, ipix)
    [z, phi] = pix2zphi(nside, order, ipix)
    [theta, phi] = pix2ang(nside, order, ipix)
    ipix = vec2pix(nside, order, vec)
    ipix = zphipix(nside, order, z, phi)
    ipix = ang2pix(nside, order, theta, phi)
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
