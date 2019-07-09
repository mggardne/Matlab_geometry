function [nx,ny,nz,xc,yc,zc] = tri_norm(tri,xyz)
%TRI_NORM Computes the normals and centroids of 3-D triangles.
%
%         [Nx,Ny,Nz,Xc,Yc,Zc] = TRI_NORM(TRI,XYZ) returns the normals,
%         [Nx,Ny,Nz] and centroids [Xc,Yc,Zc] of triangles defined by a
%         three (3) column connectivity matrix, TRI, and X, Y and Z
%         coordinates in a three (3) column matrix, XYZ.
%
%         NOTE:  1.  Must have more than one triangle.
%
%         21-Sep-2010 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error(' *** ERROR in TRI_NORM:  Two input arguments are required!');
end
%
% Check Inputs
%
[nrow1,ncol1] = size(tri);
ncol2 = size(xyz,2);
%
if (ncol1~=3)&(ncol2~=3)
  error([' *** ERROR in TRI_NORM:  Input matrices must have three', ...
        ' (3) columns!']);
end
%
% Calculate Triangle Normals
%
x = xyz(:,1)';
y = xyz(:,2)';
z = xyz(:,3)';
%
xt = x(tri);
yt = y(tri);
zt = z(tri);
%
a = [xt(:,3)-xt(:,2) yt(:,3)-yt(:,2) zt(:,3)-zt(:,2)];     % One edge
b = [xt(:,1)-xt(:,2) yt(:,1)-yt(:,2) zt(:,1)-zt(:,2)];     % Second edge
n = cross(a,b,2);                      % Normal
n = n./repmat(sqrt(sum(n.^2,2)),1,3);  % Unit normals
nx = n(:,1);
ny = n(:,2);
nz = n(:,3);
%
% Calculate Triangle Centers
%
xc = mean(xt,2);
yc = mean(yt,2);
zc = mean(zt,2);
%
return