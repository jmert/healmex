function nside = npix2nside(npix)
% nside = npix2nside(npix)
%
% INPUTS
%   npix      The number of pixels in a HEALPix map.
%
% OUTPUTS
%   nside     The corresponding HEALPix resolution factor, Nside.
%
% EXAMPLE
%   nside = healmex.npix2nside(size(map, 1));

  nside = sqrt(npix / 12);
end
