function a = rotxyz(b)
%ROTXYZ Computes the xyz-convention rotation transformation matrix.
%       ROTXYZ(B) returns the rotation matrix based on three Euler
%       angles in B.  B must be of length 3.  This convention is
%       also know as Tait-Bryan angles or the 321 sequence.
%
%       Reference:
%            Goldstein H: Classical Mechanics. 2nd Ed. Addison-Wesley
%            Publishing Co., Reading, Mass., 1980, pp 143-8, 608-10
%
%       31-Aug-94
%

%#######################################################################
%
% Check Length of Euler Angles
%
if (prod(size(b))==3)
%
% Z-Axis Rotation
%
  c = cos(b(3));
  s = sin(b(3));
  D = [c s 0; -s c 0; 0 0 1];
%
% Y-Axis Rotation
%
  c = cos(b(2));
  s = sin(b(2));
  C = [c 0 -s; 0 1 0; s 0 c];
%
% X-Axis Rotation
%
  c = cos(b(1));
  s = sin(b(1));
  B = [1 0 0; 0 c s; 0 -s c];
%
% Return Rotation Transformation Matrix
%
  a = (B*C*D)';
else
%
% Error
%
  error(' *** Error in ROTXYZ:  Input vector must be of length 3.');
%
end
%
return
