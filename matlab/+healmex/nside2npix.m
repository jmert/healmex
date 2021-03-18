function npix = nside2npix(nside)
% npix = nside2npix(nside)
%
% INPUTS
%   nside     HEALPix Nside resolution factor.
%
% OUTPUTS
%   npix      The number of pixels in a HEALPix map at the given resolution.
%
% EXAMPLE
%   npix = healmex.nside2npix(512);

  npix = 12 * nside * nside;
end
