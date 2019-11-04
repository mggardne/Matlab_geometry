function [xi,yi,zi] = gridtri(x,y,z,xi,yi,tri,method)
%GRIDTRI Data gridding and surface fitting.
%   ZI = GRIDTRI(X,Y,Z,XI,YI,TRI) fits a surface of the form Z = F(X,Y)
%   to the data in the (usually) nonuniformly-spaced vectors (X,Y,Z).
%   GRIDTRI interpolates this surface at the points specified by (XI,YI)
%   to produce ZI.  The surface is based on a triangulation of the data.
%   The surface always goes through the data points.  XI and YI are
%   usually a uniform grid (as produced by MESHGRID) and is where
%   GRIDTRI gets its name.
%
%   XI can be a row vector, in which case it specifies a matrix with
%   constant columns. Similarly, YI can be a column vector and it specifies
%   a matrix with constant rows.
%
%   [XI,YI,ZI] = GRIDTRI(X,Y,Z,XI,YI,TRI) also returns the XI and YI
%   formed this way (the results of [XI,YI] = MESHGRID(XI,YI)).
%
%   [...] = GRIDTRI(X,Y,Z,XI,YI,TRI,METHOD) where METHOD is one of
%       'linear'    - Triangle-based linear interpolation (default)
%       'cubic'     - Triangle-based cubic interpolation
%       'nearest'   - Nearest neighbor interpolation
%   defines the type of surface fit to the data. The 'cubic' method
%   produces smooth surfaces while 'linear' and 'nearest' have
%   discontinuities in the first and zero-th derivative respectively.
%   All the methods are based on a triangulation of the data.
%
%   If METHOD is [], then the default 'linear' method will be used.
%
%   Example:
%      rand('seed',0)
%      x = rand(100,1)*4-2; y = rand(100,1)*4-2; z = x.*exp(-x.^2-y.^2);
%      ti = -2:.25:2; 
%      [xi,yi] = meshgrid(ti,ti);
%      zi = gridtri(x,y,z,xi,yi,tri);
%      mesh(xi,yi,zi), hold on, plot3(x,y,z,'o'), hold off
%
%   See also GRIDDATA, GRIDDATAN, DELAUNAY, INTERP2, MESHGRID, DELAUNAYN.
%
%   * 14-Feb-2011 * Mack Gardner-Morse * Modification of griddata.
%
%   NOTE:  Function no longer checks for duplicate input coordinates.
%

%   Original copyright of griddata:
%   Copyright 1984-2004 The MathWorks, Inc. 
%   $Revision: 5.33.4.5 $  $Date: 2004/12/06 16:35:45 $

error(nargchk(6,7,nargin))

[msg,x,y,z,xi,yi] = xyzchk(x,y,z,xi,yi);
if ~isempty(msg), error(msg); end

if ( nargin < 7 || isempty(method) ),  method = 'linear'; end
if ~ischar(method), 
  error('MATLAB:gridtri:InvalidMethod',...
        'METHOD must be one of ''linear'',''cubic'', or ''nearest''.');
end

switch lower(method),
  case 'linear'
    zi = linear(x,y,z,xi,yi,tri);
  case 'cubic'
    zi = cubic(x,y,z,xi,yi,tri);
  case 'nearest'
    zi = nearest(x,y,z,xi,yi,tri);
  otherwise
    error('MATLAB:gridtri:UnknownMethod', 'Unknown method.');
end
  
if nargout<=1, xi = zi; end


%------------------------------------------------------------
function zi = linear(x,y,z,xi,yi,tri)
%LINEAR Triangle-based linear interpolation

%   Reference: David F. Watson, "Contouring: A guide
%   to the analysis and display of spacial data", Pergamon, 1994.

siz = size(xi);
xi = xi(:); yi = yi(:); % Treat these as columns
x = x(:); y = y(:); % Treat these as columns
    
if isempty(tri),
  warning('MATLAB:gridtri:NoTriangulation','Triangulation not found.');
  zi = repmat(NaN,size(xi));
  return
end

% Find the nearest triangle (t)
t = tsearch(x,y,tri,xi,yi);

% Only keep the relevant triangles.
out = find(isnan(t));
if ~isempty(out), t(out) = ones(size(out)); end
tri = tri(t,:);

% Compute Barycentric coordinates (w).  P. 78 in Watson.
del = (x(tri(:,2))-x(tri(:,1))) .* (y(tri(:,3))-y(tri(:,1))) - ...
      (x(tri(:,3))-x(tri(:,1))) .* (y(tri(:,2))-y(tri(:,1)));
w(:,3) = ((x(tri(:,1))-xi).*(y(tri(:,2))-yi) - ...
          (x(tri(:,2))-xi).*(y(tri(:,1))-yi)) ./ del;
w(:,2) = ((x(tri(:,3))-xi).*(y(tri(:,1))-yi) - ...
          (x(tri(:,1))-xi).*(y(tri(:,3))-yi)) ./ del;
w(:,1) = ((x(tri(:,2))-xi).*(y(tri(:,3))-yi) - ...
          (x(tri(:,3))-xi).*(y(tri(:,2))-yi)) ./ del;
w(out,:) = zeros(length(out),3);

z = z(:).'; % Treat z as a row so that code below involving
            % z(tri) works even when tri is 1-by-3.
zi = sum(z(tri) .* w,2);

zi = reshape(zi,siz);

if ~isempty(out), zi(out) = NaN; end
%------------------------------------------------------------

%------------------------------------------------------------
function zi = cubic(x,y,z,xi,yi,tri)
%TRIANGLE Triangle-based cubic interpolation

%   Reference: T. Y. Yang, "Finite Element Structural Analysis",
%   Prentice Hall, 1986.  pp. 446-449.
%
%   Reference: David F. Watson, "Contouring: A guide
%   to the analysis and display of spacial data", Pergamon, 1994.

if isempty(tri), 
  warning('MATLAB:gridtri:NoTriangulation','Triangulation not found.');
  zi = repmat(NaN,size(xi));
  return
end

% Find the nearest triangle (t)
t = tsearch(x,y,tri,xi,yi);

zi = cubicmx(x,y,z,xi,yi,tri,t);
%------------------------------------------------------------

%------------------------------------------------------------
function zi = nearest(x,y,z,xi,yi,tri)
%NEAREST Triangle-based nearest neightbor interpolation

%   Reference: David F. Watson, "Contouring: A guide
%   to the analysis and display of spacial data", Pergamon, 1994.

siz = size(xi);
xi = xi(:); yi = yi(:); % Treat these a columns
x = x(:); y = y(:); z = z(:); % Treat these as columns

if isempty(tri), 
  warning('MATLAB:gridtri:NoTriangulation','Triangulation not found.');
  zi = repmat(NaN,size(xi));
  return
end

% Find the nearest vertex
k = dsearch(x,y,tri,xi,yi);
zi = k;
d = find(isfinite(k));
zi(d) = z(k(d));
zi = reshape(zi,siz);
%----------------------------------------------------------