function [pxyz,norm_vec,coeff,v,score,pexp,res,sse] = plane_fit(x,y,z);
%PLANE_FIT Fits a plane to X, Y and Z point coordinate data using SVD.
%
%          [PXYZ,NORM_VEC,COEFF] = PLANE_FIT(X,Y,Z) given the X, Y
%          and Z coordinates of a set of points, calculates a least
%          squares fit of a plane. A point on the plane, PXYZ, the
%          normal vector, NORM_VEC and the coefficients of the
%          algebraic equation for a plane, COEFF, are returned.
%          The algebraic form is:
%          coeff(1)*x + coeff(2)*y + coeff(3)*z + coeff(4) = 0.
%
%          [PXYZ,NORM_VEC,COEFF,V,SCORE,PEXP,RES,SSE] = PLANE_FIT(X,Y,Z)
%          returns the rotation matrix, V, PCA scores, SCORE, percent
%          of variance explained by the three orthogonal directions of
%          the plane, PEXP, the residuals (difference between the data
%          and fitted plane), RES, and the sum of squared errors, SSE.
%
%          NOTES:  1.  Must have at least three (3) points.
%
%                  2.  The SVD is used to do a principal component
%                  analysis (PCA) to do an orthogonal regression (total
%                  least squares) fit of the plane.
%
%                  3.  Based on the demonstration of orthogonal
%                  regression using Matlab Statistics Toolbox.  See:
%                  http://www.mathworks.com/products/statistics/
%                  demos.html?file=/products/demos/shipping/stats/
%                  orthoregdemo.html (dead link).  See:
%                  http://www.mathworks.com/help/stats/examples/
%                  fitting-an-orthogonal-regression-using-principal-
%                  components-analysis.html?prodcode=ST&language=en
%
%          22-June-2010 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error([' *** ERROR in PLANE_FIT:  The X, Y and Z coordinates of ', ...
         'the points to be fit are required as inputs!']);
end
%
% Get Data Matrix and Check Number of Points
%
xyz = [x(:) y(:) z(:)];
npts = size(xyz,1);
if npts<3
  error(' *** ERROR in PLANE_FIT:  Not enough data points!');
end
%
% Center Data
%
pxyz = mean(xyz);       % Point on the fitted plane
xyz = xyz-repmat(pxyz,npts,1);         % Center data
%    
% Fit Plane
%
[u,s,v] = svd(xyz);     % Number of datapoints x 3 plane parameters
%
norm_vec = v(:,3);      % Normal vector
%
% Solve for the Coefficients
%
coeff = [norm_vec; -pxyz*norm_vec];
%
% Additional Outputs
%
if nargout>3
  score = u*s;
%
  rts = diag(s);
  pexp = 100*rts./sum(rts);
%
  res = xyz-score(:,1:2)*v(:,1:2)';
%
  err = xyz*norm_vec;
  sse = sum(err.^2);
end
%
return