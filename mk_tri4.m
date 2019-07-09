function [tri,nt,slx] = mk_tri4(dat,tol,iplt);
%MK_TRI4 Makes a triangular mesh by using the ordered slice data from
%        the digitized MRI.
%
%        [TRI,NT] = MK_TRI4(DAT) given a cell array containing three (3)
%        columns matrices with slice coordinate point data, DAT, returns
%        the three (3) column triangle connectivity matrix, TRI.  The
%        number of returned triangles, NT, may also be returned.
%
%        NOTES:  1.  Each slice coordinate data matrix must correspond
%                to one index into the cell array DAT.
%
%                2.  The coordinates should be ordered in the same
%                direction in every slice.  The dot product of the
%                directions of adjacent slices are used to check the
%                ordering direction and the ordering direction is
%                reversed if the dot product is negative.
%
%                3.  The arclength along each slice is used to determine
%                the triangulation.  See mk_tris.m for a triangulation
%                based on each slice starting at the same in slice
%                point.  See mk_triq.m for a quicker, but different
%                triangulation based on the number of points in each
%                slice.
%
%                4.  The M-file plane_fit.m must be in the current
%                path or directory.
%
%        22-Oct-2010 * Mack Gardner-Morse
%
%        16-Mar-2011 * Mack Gardner-Morse * Updated to use plane_fit for
%                                           slice separation.
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  iplt = false;
end
%
if iplt
  hf = figure;
  orient tall;
end
%
if (nargin<2)
  tol = 0.1;            % Tolerance on dot product (projection of directions)
end
%
if (nargin<1)
  error(' *** ERROR in MK_TRI4:  No input data!');
end
%
% Get Arc Lengths
%
dat = dat(:);
nslice = size(dat,1);
slen = cell(nslice,1);
npts = zeros(nslice,1);
rpt1 = zeros(nslice,3);
vec1 = zeros(nslice,3);
cntr = zeros(nslice,3);
nvec = zeros(3,nslice);
irev = false;
for k = 1:nslice
   xyz = dat{k};
   [cntr(k,:),nvec(:,k)] = plane_fit(xyz(:,1),xyz(:,2),xyz(:,3));
   vec = xyz(end,:)-xyz(1,:);
   vec1(k,:) = vec./norm(vec);
%
% Check for Slices with a Reverse Digitization
%
   if k>1
     dotp = vec1(k-1,:)*vec1(k,:)';
     if dotp<tol
       irev = true;
       xyz = flipud(xyz);
       vec = xyz(end,:)-xyz(1,:);
       vec1(k,:) = vec./norm(vec);
       dotp2 = vec1(k-1,:)*vec1(k,:)';
       if dotp2<dotp    % Revert back to original ordering
         warning([' *** WARNING in mk_tri4:  Ordering of points', ...
                  ' in the slices may not be in the same direction!']);
         irev = false;
         xyz = flipud(xyz);
         vec = xyz(end,:)-xyz(1,:);
         vec1(k,:) = vec./norm(vec);
       end
     else
       irev = false;
     end
   end
   rpt1(k,:) = xyz(1,:);
   npts(k) = size(xyz,1);
   dd = diff(xyz);
   dlen = sqrt(sum(dd.*dd,2));
   if irev
     slen{k} = flipud([0; cumsum(dlen)]);
   else
     slen{k} = [0; cumsum(dlen)];
   end
end
%
n = [0; cumsum(npts)];
tri = [];
slx = zeros(nslice-1,1);
for k = 2:nslice
%
% Slice Separations and Offsets
%
   ds = rpt1(k,:)-rpt1(k-1,:);
   slx(k-1) = (cntr(k,:)-cntr(k-1,:))*nvec(:,k-1);
   offst = ds*vec1(k-1,:)';
%
% Delaunay Triangulation
%
   xt = [zeros(npts(k-1),1); slx(k-1)*ones(npts(k),1)];
   yt = [slen{k-1}-offst; slen{k}];
 %  xt = [zeros(npts(k-1),1); ones(npts(k),1)];
 %  yt = [(0:1/(npts(k-1)-1):1)'; (0:1/(npts(k)-1):1)'];
   tril = delaunay(xt,yt);
   nid = n(k-1)+1:n(k+1);
   tri = [tri; nid(tril)];
   if iplt
     ntril = size(tril,1);
     cla;
     plot(xt,yt,'k.');
     hold on;
     trimesh(tril,xt,yt);
     text(xt,yt,int2str((1:length(xt))'),'Color','b','FontSize',12);
     text(mean(xt(tril),2),mean(yt(tril),2),int2str((1:ntril)'), ...
          'Color','r','FontSize',12);
     pause;
   end
%
end
%
nt = size(tri,1);
%
return