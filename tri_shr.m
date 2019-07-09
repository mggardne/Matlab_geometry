function hp = tri_shr(tri,xyz,shr);
%TRI_SHR  Plots 3-D triangles slightly shrunken for visualizing
%         individual triangles.
%
%         TRI_SHR(TRI,XYZ) Plots the triangles defined by the three
%         column triangle connectivity matrix, TRI, and the X, Y and Z
%         coordinates in a three column matrix, XYZ.  
%
%         HP = TRI_SHR(TRI,XYZ) Returns a patch handle for the plotted
%         triangles.
%
%         TRI_SHR(TRI,XYZ,SHR) Optional variable, SHR, controls the
%         amount the triangles are shrunk.  The variable should be
%         greater than 0 (plots just the centers of the triangles) and
%         1 (triangles not shrunk at all).  By default SHR = 0.75.
%         Generally, SHR should be between 0.5 and 0.85.
%
%         NOTES:  1.  Color of triangular patches based on the Z
%                 coordinate (height).
%
%         16-May-2012 * Mack Gardner-Morse
%

%#######################################################################
%
% Check the Inputs
%
if (nargin<2)
  error(' *** ERROR in TRI_SHR:  Not enough input arguments!');
end
%
if (nargin<3)
  shr = 0.75;
end
%
if shr<=0|shr>1
  shr = 0.75;
end
%
% Number of Triangles
%
nt = size(tri,1);
%
% Get Shrunken Coordinates
%
xp = reshape(xyz(tri,1),nt,3)';
yp = reshape(xyz(tri,2),nt,3)';
zp = reshape(xyz(tri,3),nt,3)';
xp = repmat(mean(xp),3,1)+shr*(xp-repmat(mean(xp),3,1));
yp = repmat(mean(yp),3,1)+shr*(yp-repmat(mean(yp),3,1));
zp = repmat(mean(zp),3,1)+shr*(zp-repmat(mean(zp),3,1));
%
% Plot Triangles
%
hp = patch(xp,yp,zp,zp);
% view(3);
% axis equal;
%
return