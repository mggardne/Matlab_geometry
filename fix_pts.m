function [xyz,npts] = fix_pts(xyz,tol,iflag);
%FIX_PTS Removes duplicate points from a matrix of 3-D data.
%
%        [XYZ,NPTS] = FIX_PTS(XYZ,TOL) given a three (3) columns matrix
%        with coordinate point data, XYZ, returns the matrix without
%        any duplicate points within a distance, TOL.  The number of
%        returned points, NPTS, may also be returned.
%
%        [XYZ,NPTS] = FIX_PTS(XYZ,TOL,IFLAG) if IFLAG is true (or not
%        equal to zero (0)), prints a message to the screen if any
%        duplicates are found.
%
%        NOTES:  None.
%
%        25-June-2010 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  iflag = false;
end
%
if (nargin<2)
  error(' *** ERROR in FIX_PTS:  Must have two inputs!');
end
%
[npts,ncol] = size(xyz);
if ncol~=3
  error(' *** ERROR in FIX_PTS:  Data must have three columns!');
end
%
% Find and Remove Duplicates
%
tol2 = tol.*tol;
npts = size(xyz,1);
for k = 1:npts-1;
   l = k+1:npts;
   dxyz = xyz(l,:)-repmat(xyz(k,:),npts-k,1);
   dist2 = sum(dxyz.*dxyz,2);
   idx = find(dist2<tol2);
   if ~isempty(idx)
     if iflag
       fprintf(1,'\n *** Duplicate Points Found!\n');
%        error('\n *** Duplicate Points Found!\n');
       iflag = false;
     end
     idx = [k; l(idx)'];
     ndup = size(idx,1);
     xyz(k,:) = mean(xyz(idx,:));
     xyz(idx(2:ndup),:) = repmat(NaN,ndup-1,3);
   end
end
%
xyz = xyz(~isnan(xyz(:,1)),:);
npts = size(xyz,1);
%
return