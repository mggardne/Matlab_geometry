function [davg,dstd,dmin,dmax,daavg,dif,xyz1,tri1,xyz2,tri2] = ...
         dig_zdif3(dat1,dat2,iplt,hf);
%DIG_ZDIF3 Gets the Z difference between two MRI digitizations of a
%          surface using a plane fit and 2D interpolation except at the
%          boundary (edge) points of the surface.
%
%         [DAVG,DSTD,DMIN,DMAX] = DIG_ZDIF3(DAT1,DAT2) given cell arrays
%         of digitized MRI surface data, DAT1 and DAT2, fits a plane,
%         interpolates the Z coordinates of DAT1 to the X, Y coordinates
%         of DAT2 and takes the differences in the Z coordinates.  The
%         mean difference, DAVG, the standard deviation of the
%         differences, DSTD, minimum difference, DMIN, (maximum negative
%         difference), and maximum difference, DMAX, are returned.
%
%         [DAVG,DSTD,DMIN,DMAX,DAAVG] = DIG_ZDIF3(DAT1,DAT2) returns
%         the mean absolute difference, DAAVG.
%
%         [DAVG,DSTD,DMIN,DMAX,DAAVG,DIF] = DIG_ZDIF3(DAT1,DAT2) returns
%         a column vector of Z differences, DIF, between the two
%         digitizations.  The differences will contain NaNs where DAT2
%         XY data is outside of the DAT1 XY data and along the boundary
%         of the DAT1 data.
%
%         [DAVG,DSTD,DMIN,DMAX] = DIG_ZDIF3(DAT1,DAT2,IPLT) if IPLT is
%         true, the surfaces and differences are plotted in a new
%         figure.  The DAT2 surface is shown higher with red bars and
%         the DAT2 surface is shown lower with blue bars.  Maximum
%         differences are shown in wider and brighter red and blue bars.
%
%         [DAVG,DSTD,DMIN,DMAX] = DIG_ZDIF3(DAT1,DAT2,IPLT,HF) the
%         the surfaces and differences are plotted in the figure with
%         figure handle HF.
%
%         NOTES:  1.  The surfaces are plotted in the transformed best-
%                 fit plane and not in the input coordinate system.
%
%                 2.  Works best on surfaces that are uniquely defined
%                 in the XY plane.
%
%                 3.  The M-files in_tri2d.m, meshbnd4.m, mk_tri4.m,
%                 plane_fit.m and tri_fix2.m must be in the current
%                 path or directory.
%
%         09-Jul-2019 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)||isempty(iplt)
  iplt = false;
end
%
if (nargin<4)||isempty(hf)
  if iplt
    hf = figure;
  end
end
%
if (nargin<2)
  error(' *** ERROR in DIG_ZDIF3:  Not enough input data!');
end
%
% Check Inputs
%
if ~iscell(dat1)||~iscell(dat2)
  error(' *** ERROR in DIG_ZDIF3:  Inputs must be cell arrays!');
end
%
dat1 = dat1(:);
dat2 = dat2(:);
%
% Get Surface Triangulations and Coordinates
%
xyz1 = cell2mat(dat1);  % First digitization
tri1 = tri_fix2(mk_tri4(dat1),xyz1);   % First digitization
%
xyz2 = cell2mat(dat2);  % Second digitization
tri2 = tri_fix2(mk_tri4(dat2),xyz2);   % Second digitization
%
% Get Best Fit Plane
%
xyz3 = [xyz1; xyz2];
[pxyz,nvec,coeff,r] = plane_fit(xyz3(:,1),xyz3(:,2),xyz3(:,3));
%
n1 = size(xyz1,1);
xyzo1 = repmat(pxyz,n1,1);
xyz1 = (xyz1-xyzo1)*r+xyzo1;
%
n2 = size(xyz2,1);
xyzo2 = repmat(pxyz,n2,1);
xyz2 = (xyz2-xyzo2)*r+xyzo2;
%
% Make Sure Normal Vector is in Postive Direction
%
if nvec(3)<0
  xyz1(:,3) = -xyz1(:,3);
  xyz2(:,3) = -xyz2(:,3);
end
%
% Get Boundary Nodes
%
bid = meshbnd4(tri2);
%
% Check Consistency of Digitization
%
tgrid = scatteredInterpolant(xyz1(:,[1 2]),xyz1(:,3),'linear','none');
z2 = tgrid(xyz2(:,[1 2]));
it = in_tri2d(tri1,xyz1(:,[1 2]),xyz2(:,[1 2]));
z2(~it) = NaN;
%
dif = z2-xyz2(:,3);
%
dif(bid) = NaN;         % Set differences on boundary to NaNs
%
[dmin,idmn] = min(dif);
[dmax,idmx] = max(dif);
%
idnnan = ~isnan(dif);
dif2 = dif(idnnan);     % Remove NaNs
dstd = std(dif2);
davg = mean(dif2);
daavg = mean(abs(dif2));
%
if iplt
  figure(hf);
  orient landscape;
  hm = trimesh(tri1,xyz1(:,1),xyz1(:,2),xyz1(:,3),'EdgeColor', ...
               [0.75 0.75 0.75],'FaceColor','none','LineWidth',0.5);
%  set(hm,'EdgeAlpha',0.5);
  axis equal;
  hold on;
%
  idn = find(dif<0);
  idp = find(dif>0);
%
  plot3([xyz2(idn,1) xyz2(idn,1)]',[xyz2(idn,2) ...
        xyz2(idn,2)]',[z2(idn) xyz2(idn,3)]','r-', ...
        'Color',[0.5 0 0],'LineWidth',1);
%  
  plot3([xyz2(idp,1) xyz2(idp,1)]',[xyz2(idp,2) ...
        xyz2(idp,2)]',[z2(idp) xyz2(idp,3)]','b-', ...
        'Color',[0 0 0.5],'LineWidth',1);
%  
  plot3([xyz2(idmn,1); xyz2(idmn,1)],[xyz2(idmn,2); ...
        xyz2(idmn,2)],[z2(idmn); xyz2(idmn,3)],'r-', ...
        'LineWidth',2);
%
  plot3([xyz2(idmx,1); xyz2(idmx,1)],[xyz2(idmx,2); ...
         xyz2(idmx,2)],[z2(idmx); xyz2(idmx,3)],'b-', ...
        'LineWidth',2);
%
  axlim = axis;
  txtc = axlim([1 4 6]);
  txtc = [0.98 1.05 1.05].*(txtc<0).*txtc+ ...
         [1.02 0.95 0.95].*(txtc>0).*txtc;
  text(txtc(1),txtc(2),txtc(3),{['Mean Difference = ', ...
       num2str(davg) ' mm']; ['SD = ' num2str(dstd) ' mm']; ...
       ['Minimum Difference = ' num2str(dmin) ' mm']; ...
       ['Maximum Difference = ' num2str(dmax) ' mm']; ...
       ['Mean Absolute Difference = ' num2str(daavg) ' mm']}, ...
       'FontSize',11,'FontWeight','bold');
  xlabel('X (mm)','FontSize',12,'FontWeight','bold');
  ylabel('Y (mm)','FontSize',12,'FontWeight','bold');
  zlabel('Z (mm)','FontSize',12,'FontWeight','bold');
  title('Digitization Differences','FontSize', ...
        16,'FontWeight','bold','Interpreter','none');
end
%
return