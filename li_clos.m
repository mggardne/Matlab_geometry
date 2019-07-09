function [p1,p2,dist] = li_clos(li1,li2);
%LI_CLOS Finds the closest two points between two nonparallel 3D lines.
%
%        [P1,P2] = LI_CLOS(LI1,LI2) given a two (2) rows by three (3)
%        columns matrix giving 3D coordinates for two points on one
%        line, LI1, and another two (2) rows by three (3) columns
%        matrix giving 3D coordinates for two points on a second line,
%        LI2, returns the three (3) column row vector of the coordinates
%        for the point on the first line that is closest to the second
%        line, P1, and the the three (3) column row vector of the
%        coordinates for the point on the second line that is closest to
%        the first line, P2.  Note that the line between P1 and P2 is
%        also perpendicular to both the input lines.
%
%        [P1,P2,DIST] = LI_CLOS(LI1,LI2) returns the closest distance
%        between the two lines, DIST.
%
%        NOTES:  1.  Based on the algorithm at:
%        http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm
%
%                2.  Does not check that the two points on the lines are
%                not the same.
%
%        17-Aug-2010 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error([' *** ERROR in LI_CLOS:  At least two (2) input lines', ...
         ' are required!']);
end
%
% Check Inputs
%
[r1,c1] = size(li1);
[r2,c2] = size(li2);
%
if r1~=2|r2~=2
  error([' *** ERROR in LI_CLOS:  Two (2) points are', ...
         ' required to define each line!']);
end
%
if c1~=3|c2~=3
  error([' *** ERROR in LI_CLOS:  Three (3) coordinates are', ...
         ' required to define each point!']);
end
%
% Get Line Direction Vectors
%
r0 = li1(1,:)-li2(1,:);
u = diff(li1);
v = diff(li2);
%
% Get Dot Products
%
a = u*u';
b = u*v';
c = v*v';
d = u*r0';
e = v*r0';
%
% Solve for Line Parameters
%
den = a*c-b*b;
tol = eps^(2/3);
if abs(den)<tol
  error([' *** ERROR in LI_CLOS:  The two lines are parallel', ...
         ' or almost parallel!']);
else
  sc = (b*e-c*d)./den;
  tc = (a*e-b*d)./den;
end
%
% Solve for the Points
%
p1 = li1(1,:)+sc*u;
p2 = li2(1,:)+tc*v;
%
% Get Shortest Distance between Lines
%
if nargout>2
  dist = p2-p1;
  dist = sqrt(dist*dist');
end
%
return