function [x,y,z] = sph2cartGu(aziDeg, eleDeg, rad)
% [x,y,z] = sph2cartGu(aziDeg, eleDeg, rad) converts a 3D point from
% spherical to cartesian coordinates. Azimuth and elevation are given in
% degrees.
% The coordinate system is laid out so that it matches the Gu et al. 
% studies (Gu et al., 2006; Takahashi et al., 2007) as well as the
% FlowRaudies2012 object.
%
% x-coordinate: rightward-leftward
% y-coordinate: downward-upward
% z-coordinate: forward-backward

[x,z,y] = sph2cart(deg2rad(aziDeg), deg2rad(eleDeg), rad);
y = -y;

end