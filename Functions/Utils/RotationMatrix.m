function [rotmat] = RotationMatrix(phi,axis)
% Usage [rotmat] = RotationMatrix(phi,axis)
% 
% Gives a 3D rotation matrix
% phi  =  the angle (in rad) of rotation
% axis =  the axis around which the rotation is performed
%
% JvB, Sept 2011

    rotmat = ...
    [ cos(phi)+axis(1)^2*(1-cos(phi))                   axis(1)*axis(2)*(1-cos(phi))-axis(3)*sin(phi)       axis(1)*axis(3)*(1-cos(phi))+axis(2)*sin(phi)   ;...
      axis(2)*axis(1)*(1-cos(phi))+axis(3)*sin(phi)     cos(phi)+axis(2)^2*(1-cos(phi))                     axis(2)*axis(3)*(1-cos(phi))-axis(1)*sin(phi)   ;...
      axis(3)*axis(1)*(1-cos(phi))-axis(2)*sin(phi)     axis(3)*axis(2)*(1-cos(phi))+axis(1)*sin(phi)       cos(phi)+axis(3)^2*(1-cos(phi))     ];
end