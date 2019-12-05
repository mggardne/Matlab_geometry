function zg = gridproj(trip,xyzp,xg,yg,zdir)
%GRIDPROJ Projects a two dimensional grid along a vector perpendicular
%         to the grid onto a triangular mesh and finds the most distant
%         intersection along the vector.
%
%         ZG = GRIDPROJ(TRIP,XYZP,XG,YG,ZDIR) given a three (3) columns
%         triangle connectivity matrix, TRIP, a three (3) columns
%         coordinate point data matrix, XYZP, grid X and Y coordinates
%         in vectors, XG and YG, and a Z direction parameter, ZDIR,
%         projects the grid onto the surface defined by the mesh and
%         returns the most distant surface intersection Z coordinate in
%         the positive Z direction if ZDIR is greater than or equal to
%         zero or in the negative Z direction if ZDIR is less than zero.
%         The most distant surface intersection Z coordinates are
%         returned in vector, ZG.
%
%         NOTES:  1.  Function to find the posterior of the bony
%                 patella.
%
%                 2.  Finds mesh points (and connecting triangles)
%                 within 2.4 units of the projected lines.  See line 63.
%
%                 3.  The M-files nod2tri.m, pts2lin.m, tsect4.m and
%                 xprod.m must be in the current path or directory.
%
%         09-Dec-2015 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<4)
  error(' *** ERROR in GRIDPROJ:  Not enough inputs!');
end
%
if (nargin<5)
  zdir = 1;             % Positive direction by default
end
%
% Check Inputs
%
ncolt = size(trip,2);
ncolp = size(xyzp,2);
%
if (ncolt~=3)||(ncolp~=3)
  error([' *** ERROR in GRIDPROJ:  Mesh inputs must have', ...
         ' three (3) columns!']);
end
%
xg = xg(:);
yg = yg(:);
ngp = size(xg,1);
ngpt = size(yg,1);
%
if ngp~=ngpt
  error([' *** ERROR in GRIDPROJ:  Grid inputs must have', ...
         ' the same lengths!']);
end
%
% Analysis Parameters
%
dist = 2.4*2.4;         % Find points within two (2) scan widths of the line
%
% Get Projection Direction
%
if zdir<0
  zdir = -1;
else
  zdir = 1;
end
%
zvec = [0 0 zdir];      % Projection vector
%
% Loop Through Grid Points
%
zg = zeros(ngp,1);      % Initialize array
%
for k = 1:ngp
%
% Find Points Within a Distance of the Grid Projection Line
%
   gpt = [xg(k) yg(k) 0];
   xyzl = pts2lin(gpt,zvec,xyzp);
   d = xyzl-xyzp;
   d = sum(d.*d,2);
   idp = find(d<dist);  % Points within sqrt(dist) of projection line
%
% Check for Intersections with Mesh Surface
%
   if ~isempty(idp)
     idt = nod2tri(idp,trip);
     nt = size(idt,1);
%
% Loop through Triangles Close to Projection Line Looking for
% Intersections
%
     xyzi = NaN(nt,3);
%
     for l = 1:nt
        it = trip(idt(l),:)';          % Nodes in triangle
        v1 = xyzp(it(1),:);
        v2 = xyzp(it(2),:);
        v3 = xyzp(it(3),:);
        [xyzit,il] = tsect4(v1,v2,v3,gpt,zvec);
        if il
          xyzi(l,:) = xyzit';          % Intersection point
        end
     end
%
% Find Most Distant Intersection Point
%
     idd = find(~isnan(xyzi(:,3)));
     if ~isempty(idd)
       zg(k) = zdir*max(zdir*xyzi(idd,3));
     else
       zg(k) = NaN;
     end
   else
     zg(k) = NaN;
   end
end
%
return