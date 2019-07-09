function [area,cg] = tri_area(x,y,z,tric);
%TRI_AREA Finds the area of 3-D triangles.
%
%         AREA = TRI_AREA(X,Y,Z,TRIC) returns the areas of the
%         triangles defined by the X, Y and Z coordinates and the
%         three (3) column triangle connectivity matrix TRIC.
%
%         [AREA,CG] = TRI_AREA(X,Y,Z,TRIC) returns the area centers of
%         the triangles defined by the X, Y and Z coordinates and the
%         three (3) column triangle connectivity matrix TRIC in the
%         three (3) column matrix CG with the X, Y and Z coordinates of
%         the area centers in the columns
%
%         NOTES:  None.
%
%         13-Dec-04 * Mack Gardner-Morse
%
%         01-Jul-2010 * Mack Gardner-Morse * Include area centers.
%
%         17-Apr-2014 * Mack Gardner-Morse * Made coordinates row
%         vectors so that the function works with just one (1) input
%         triangle.
%

%#######################################################################
%
% Check Input Parameters
%
if (nargin<4)
  error(' *** ERROR in TRI_AREA:  Not enough input arguments.');
end
%
% Make Sure Coordinates are Row Vectors
%
x = x(:)';
y = y(:)';
z = z(:)';
%
% Find Areas
%
xe = x(tric);
ye = y(tric);
ze = z(tric);
%
v1 = [xe(:,2)-xe(:,1) ye(:,2)-ye(:,1) ze(:,2)-ze(:,1)]';
v2 = [xe(:,3)-xe(:,1) ye(:,3)-ye(:,1) ze(:,3)-ze(:,1)]';
%
c = cross(v1,v2);
%
area = sqrt(sum(c.*c))';
area = area/2;
%
% Find Centers
%
if nargout>1
  cg = [mean(xe,2) mean(ye,2) mean(ze,2)];
end
%
return