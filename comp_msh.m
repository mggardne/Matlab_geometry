function [xyzn,trin] = comp_msh(xyzo,trio);
%COMP_MSH Reduces the number of nodal coordinates to the number of
%         nodes in the mesh and renumbers the mesh connectivity to
%         match the compressed number of nodes.
%
%         [XYZN,TRIN] = COMP_MSH(XYZO,TRIO) given the X, Y and Z
%         coordinates in the columns of matrix XYZO and a mesh
%         connectivity matrix TRIO, reduces the number of nodal
%         coordinates to the number of nodes in the mesh which is
%         returned in XYZN.  Also, renumbers the mesh connectivity to
%         match the compressed number of nodes which is returned in
%         TRIN.
%
%         NOTES:  1.  For reducing the sizes of submesh coordinates.
%
%         26-Apr-2013 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error(' *** ERROR in COMP_MSH:  COMP_MSH requires two inputs!');
end
%
% Reduce the Number of Nodes and Renumber Connectivity
%
nnew = unique(trio(:));
xyzn = xyzo(nnew,:);
%
nnod = size(xyzn,1);
trin = trio;
for k = 1:nnod
   trin(trio==nnew(k)) = k;
end
%
return