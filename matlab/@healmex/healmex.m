classdef healmex < matlab.mixin.CustomDisplay
  properties (Constant, Access = private)
    id_heartbeat        = int64(-1)
    id_pix2ang          = int64(3)
  end
  methods (Static)
    [theta,phi] = pix2ang(nside, ipix);
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
