function [aziDeg,eleDeg,rad] = cart2sphGu(x, y, z)
% [aziDeg,eleDeg,rad] = cart2sphGu(x, y, z) converts a 3D point from
% cartesian to spherical coordinates. Azimuth and elevation are given in
% degrees.
% The coordinate system is laid out so that it matches the Gu et al. 
% studies (Gu et al., 2006; Takahashi et al., 2007) as well as the
% FlowRaudies2012 object.
%
% x-coordinate: rightward-leftward
% y-coordinate: downward-upward
% z-coordinate: forward-backward

[azi,ele,rad] = cart2sph(x, z, -y);

% azimuth needs to be corrected from [-180,180] to [0,360]
azi(azi<0) = azi(azi<0) + 2*pi;
aziDeg = rad2deg(azi);

% elevation is fine between [-90,90]
eleDeg = rad2deg(ele);

end