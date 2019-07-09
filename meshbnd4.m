function out = meshbnd4(tri);
%MESHBND4 Finds and orders a list of boundary nodes from lowest to 
%         highest node number.  See MESHBND3.M for a similar method
%         where nodes are ordered to form a closed connected boundary
%         around the outside edges of a triangular mesh.
%
%         OUT = MESHBND4(TRI) given a three (3) column triangle
%         connectivity matrix, TRI, returns an ordered list of boundary
%         nodes, OUT, from lowest to highest node number.
%
%         NOTES:  1.  The node list starts with the lowest node ID and
%                 ends with the highest node ID of the boundary nodes.
%
%                 2.  See MESHBND3.M and MESHBND2.M for similar methods
%                 where nodes are ordered to form a closed connected
%                 boundary around the outside edges of a triangular
%                 mesh.
%
%         01-Aug-2013 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<1
  error([' *** ERROR in MESHBND4:  A triangle connectivity matrix' ...
         ' is required as an input!']);
end
%
[nr,nc] = size(tri);
if nc~=3
  error([' *** ERROR in MESHBND4:  Triangle connectivity matrix' ...
         ' must have three columns!']);
end
%
% Get Sides of the Triangles
%
sides = [tri(:,1:2); tri(:,2:3); tri(:,[3 1])];
sides = sort(sides,2);
ns = size(sides,1);
%
% Find Duplicate Edges (Duplicate Edge == Interior Edge)
%
sides = sortrows(sides);
ds = diff(sides);
ds = sum(abs(ds),2);
idup = find(ds==0);
%
% Delete Duplicate Edges
%
if ~isempty(idup)
  idup = [idup; idup+1];
  idx = true(ns,1);
  idx(idup) = false;
  sides = sides(idx,:);
end
%
% Get Unique and Sorted Node IDs
%
out = unique(sides(:));
%
return