function tri = tri_fix2(tri,xyz);
%TRI_FIX2 Trys to improve the aspect ratios of triangles in a
%         three-dimensional (3-D) triangular mesh.
%
%         TRI = TRI_FIX2(TRI,XYZ) given a three column triangular
%         connectivity matrix TRI and the X, Y and Z coordinates in the
%         columns of matrix XYZ, checks the sum of the opposite angles
%         of the triangles sharing a common edge and flips the edge if
%         the sum of the opposite angles is greater than pi.
%
%         NOTES:  1.  Only operates locally (pairs of triangles) to
%                 improve the aspect ratio.  Other more global
%                 functions or functions that change the number of
%                 triangles may work better to improve the mesh.
%
%                 2.  The coordinates of the vertices of the triangles
%                 must be a matrix with the vertices in the rows and
%                 the X, Y and Z coordinates in the columns.
%
%                 3.  The M-file tri_norm.m must be in the current path
%                 or directory.
%
%         27-Aug-2013 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error(' *** ERROR in TRI_FIX2:  TRI_FIX2 requires two inputs!');
end
%
% Check Inputs
%
[nt,nc1] = size(tri);
[nc,nc2] = size(xyz);
if nc1~=3||nc2~=3
  error([' *** ERROR in TRI_FIX2:  Triangle connectivity matrix' ...
         ' and coordinate matrix must have three columns!']);
end
%
% Get Edges
%
edges = [tri(:,1:2); tri(:,2:3); tri(:,[3 1])];
edges = sort(edges,2);                 % Order vertices within edges
%
ne = size(edges,1);                    % Number of edges
idt = repmat((1:nt)',3,1);             % Triangles corresponding to edges
%
% Find Duplicate or Common Edges
%
[edges,is] = sortrows(edges);          % Sort edges
idt = idt(is);                         % Triangles corresponding to sorted edges
ids = repmat(1:3,nt,1);
ids = ids(:);
ids = ids(is);                         % Triangle sides corresponding to sorted edges
%
ds = diff(edges);
ds = sum(abs(ds),2);
idup = find(ds==0);                    % Duplicate or common edges
nd = size(idup,1);                     % Number of common edges
%
idc = [idup idup+1]';
idc = idc(:);
edges = edges(idc,:);
idt = idt(idc);
ids = ids(idc);
%
% Get Scores for Opposite Angles
%
sc = zeros(nd,1);
%
for k = 1:nd
%
   l = 2*k;
   n23 = edges(l,:);                   % Common edge
%
   n1 = setdiff(tri(idt(l-1),:),n23);
   p1 = xyz(n1,:);
   v1 = xyz(n23(1),:)-p1;
   v2 = xyz(n23(2),:)-p1;
   ang1 = v1*v2';                      % Dot product
   ang1 = ang1./sqrt(v1*v1'*v2*v2');
   ang1 = acos(ang1);
%
   n2 = setdiff(tri(idt(l),:),n23);
   p2 = xyz(n2,:);
   v3 = xyz(n23(1),:)-p2;
   v4 = xyz(n23(2),:)-p2;
   ang2 = v3*v4';                      % Dot product
   ang2 = ang2./sqrt(v3*v3'*v4*v4');
   ang2 = acos(ang2);
   sc(k) = ang1+ang2-pi;               % Scores > 0 should be flipped
end
%
% Sort Scores
%
[ssc,iss] = sort(sc,1,'descend');
%
% Flip Edges and Update Score
%
iter = 0;
while ssc(1)>0
     iter = iter+1;
%
% Form New Triangles by Flipping Edge
%
     l = 2*iss(1);
     n23 = edges(l,:);                 % Common edge
     idt1 = idt(l-1);
     n1 = setdiff(tri(idt1,:),n23);
     idt2 = idt(l);
     n2 = setdiff(tri(idt2,:),n23);
     nid1 = [n1 n23(1) n2];
     nid2 = [n2 n23(2) n1];
%
% Check Node Ordering of the Triangles to Preserve Triangles Normal
% Directions
%
     [xn,yn,zn] = tri_norm(tri([idt1;idt2],:),xyz);
     nv = [mean(xn) mean(yn) mean(zn)];          % Mean normal vector
     nv = nv/norm(nv);
%
     [xn,yn,zn] = tri_norm(nid1,xyz);
     nvn = [xn,yn,zn];
     ang = acos(nvn*nv');
     if ang>pi/2
       nid1 = nid1([1 3 2]);           % Reverse ordering to flip normal
     end
%
     [xn,yn,zn] = tri_norm(nid2,xyz);
     nvn = [xn,yn,zn];
     ang = acos(nvn*nv');
     if ang>pi/2
       nid2 = nid2([1 3 2]);           % Reverse ordering to flip normal
     end
%
% Update Triangle Connectivity Matrix, Edges, Triangle Index and Side Index
%
     tri(idt1,:) = nid1;
     tri(idt2,:) = nid2;
%
     ne1 = [tri(idt1,1:2); tri(idt1,2:3); tri(idt1,[3 1])];
     ne2 = [tri(idt2,1:2); tri(idt2,2:3); tri(idt2,[3 1])];
%
     ne1 = sort(ne1,2);                % Order vertices within edges
     ne2 = sort(ne2,2);                % Order vertices within edges
%
     [s,i1,i2] = intersect(ne1,ne2,'rows');
     edges(l-1:l,:) = repmat(s,2,1);
     idt(l-1:l) = [idt1;idt2];
     ids(l-1:l) = [i1;i2];
%
     [s1,i1] = setdiff(ne1,ne2,'rows');
     [s2,i2] = setdiff(ne2,ne1,'rows');
     si = [i1; i2];
     [s,i3,i4] = intersect(edges,[s1;s2],'rows');
     i3 = i3+rem(i3,2);
     idx = [i3-1 i3]';
     idx = idx(:);
     idx2 = find(idt(idx)==idt1|idt(idx)==idt2);
     idx2 = idx(idx2);
     idts = [idt1; idt1; idt2; idt2];
     idts = idts(i4);
% iter
% if iter>=28; keyboard; end
     idt(idx2) = idts;
     si = si(i4);
     ids(idx2) = si;
%
% Update Scores of Adjoining Edges
%
     sc(iss(1)) = -sc(iss(1));
     ids1 = ids(l-1);
     ids2 = ids(l);
     ide1 = find(idt==idt1&ids~=ids1);
     ide2 = find(idt==idt2&ids~=ids2);
     ide = unique([ide1; ide2]);
     ide = ide+rem(ide,2);
     ns = size(ide,1);
     if ns>0
       for k = 1:ns
          l = ide(k);
          n23 = edges(l,:);            % Common edge
%
          n1 = setdiff(tri(idt(l-1),:),n23);
          p1 = xyz(n1,:);
          v1 = xyz(n23(1),:)-p1;
          v2 = xyz(n23(2),:)-p1;
          ang1 = v1*v2';               % Dot product
          ang1 = ang1./sqrt(v1*v1'*v2*v2');
          ang1 = acos(ang1);
%
          n2 = setdiff(tri(idt(l),:),n23);
          p2 = xyz(n2,:);
          v3 = xyz(n23(1),:)-p2;
          v4 = xyz(n23(2),:)-p2;
          ang2 = v3*v4';               % Dot product
          ang2 = ang2./sqrt(v3*v3'*v4*v4');
          ang2 = acos(ang2);
          m = l/2;
          sc(m) = ang1+ang2-pi;        % Scores > 0 should be flipped
       end
     end
%
% Get New Sort Order for Scores
%
     [ssc,iss] = sort(sc,1,'descend');
end
%
return