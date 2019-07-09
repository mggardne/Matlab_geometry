function [xyzo,xyzax,aspect,widt,height] = tibia_cs8(faxial,leg,iplt)
%TIBIA_CS8 Reads Osirix regions-of-interest comma separated files to
%          calculate tibia based coordinate systems (CS).
%
%          [XYZO,XYZAX] = TIBIA_CS8(FAXIAL,LEG) given the filename of
%          an Osirix regions-of-interest comma separated of transverse
%          sections of the tibia, FAXIAL, and an integer 1 (or logical
%          true) for a right knee and 0 (or logical false) for a left
%          knee, LEG, calculates tibia coordinate system based on the
%          centroids from the most distal section of the tibia visible
%          on the MRI and the most proximal section of the tibia.  The
%          most posterior points on the tibia at the transverse section
%          through the PCL form the second axis.  The function returns 
%          the X, Y and Z center of the tibia coordinate system, XYZO,
%          and unit vectors for the X, Y and Z axes in the columns of
%          matrix XYZAX.
%
%          [XYZO,XYZAX] = TIBIA_CS8(FAXIAL,LEG,IPLT) generates plots as
%          it calculates the axes.
%
%          [XYZO,XYZAX,ASPECT,WIDT,HEIGHT] = TIBIA_CS8(FAXIAL,LEG)
%          returns the aspect ratio, ASPECT, of the width of the tibia
%          plateau, WIDT, to the height between the distal and proximal
%          outlines of the tibia plateau, HEIGHT.
%
%          NOTES:  1.  The filename must be in single row character
%                  array.
%
%                  2.  The names of the regions of interest (ROIs) must
%                  follow a naming convention.
%
%                  3.  The name of the transverse ROI of the posterior
%                  tibia must be "postaxis".
%
%                  4.  The names of the distal transverse ROIs of the
%                  tibias must be "dist_tibia" and the proximal
%                  sections must be named "prox_tibia".
%
%                  5.  The Matlab files li_clos.m, rd_roi4.m, rotxyz.m
%                  and tri_area.m must be in the current directory or
%                  path.
%
%                  6.  This is a cleaned up version of tibia_cs7.m with
%                  a loop to find improved posterior points of the
%                  tibia.
%
%          07-May-2019 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error([' *** ERROR in tibia_cs8:  At least an input', ...
         ' filename and left/right knee identifier are required!']);
end
%
if (nargin<3)
  iplt = false;
end
%
% File Name
%
nrow = size(faxial,1);
if nrow~=1
  error([' *** ERROR in TIBIA_CS8:  Filename must be in single', ...
         ' row character array!']);
end
%
% Get Current Figure
%
if iplt
  gcf;
  gca;
  hold on;
end
%
% Posterior Points Parameters
%
ptol = pi/180;          % One degree in radians
ang = zeros(2,1);       % Rotation of Y-axis based on posterior points
dang = 10*ptol;         % Initial difference in rotation angle
%
% Initialize Variables
%
xme = zeros(2,1);       % Most posterior points on PCL slice
yme = zeros(2,1);       % Most posterior points on PCL slice
zme = zeros(2,1);       % Most posterior points on PCL slice
%
% Read File
%
filenam = deblank(faxial);
%
% Get ROI Data
%
roi = rd_roi4(filenam);
nroi = size(roi,1);
%
% Check for DIST_TIBIA ROI
%
iuse = false;           % Use dist_tibia ROI and NOT anataxis ROI
%
roinams = lower(char(roi.name));
idd = strmatch('dist_tibia',roinams);
if isempty(idd)
  idd = strmatch('anataxis',roinams);
  if isempty(idd)
    error([' *** ERROR in tibia_cs8:  No distal tibia ROI found in', ...
           ' data file!']);
  end
  iuse = true;          % Use anataxis ROI and NOT missing dist_tibia ROI
end
%
% Loop through ROIs
%
for l = 1:nroi
%
   roiname = roi(l).name;
   roinaml = lower(roiname);
%
   dat = roi(l).data';
   nslice = size(dat,1);
%
% Get Distal Tibia Centroids
%
   if strcmp(roinaml,'anataxis')&iuse
%
% Get Digitized Data
%
     xyzc = zeros(nslice,3);
%
     for n = 1:nslice
        xyz = dat{n};
        npts = size(xyz,1);
%
        xyzm = mean(xyz);
        xyzp = [xyz; xyzm];
%
% Get Centroids
%
        tri = [repmat(npts+1,npts-1,1) (1:npts-1)' (2:npts)'];
        tri = [tri; [npts+1 npts 1]];
%
        [at,cgt] = tri_area(xyzp(:,1),xyzp(:,2),xyzp(:,3),tri);
        xyzc(n,:) = sum(repmat(at,1,3).*cgt)./sum(at);
        
        if iplt
          plot3(xyzc(n,1),xyzc(n,2),xyzc(n,3),'bs', ...
                'LineWidth',1,'MarkerSize',8);
        end
     end
%
% Get and Plot Most Distal Centroid
%
     [~,idz] = min(xyzc(:,3));
     xyzc1 = xyzc(idz,:);              % Most distal centroid
     if iplt
       xyz = dat{idmn};
     end
%
     if iplt
       plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k.-','MarkerSize',8, ...
             'LineWidth',1);

       plot3(xyzc1(:,1),xyzc1(:,2),xyzc1(:,3),'bs','MarkerSize',8, ...
             'LineWidth',1);
     end
%
   end
%
% Get Most Distal Tibia Centroid
%
   if strcmp(roinaml,'dist_tibia')
%
     xyz = cell2mat(dat);
     npts = size(xyz,1);
%
     xyzm = mean(xyz);
     xyzp = [xyz; xyzm];
%
% Get and Plot Centroid
%
     tri = [repmat(npts+1,npts-1,1) (1:npts-1)' (2:npts)'];
     tri = [tri; [npts+1 npts 1]];
%
     [at,cgt] = tri_area(xyzp(:,1),xyzp(:,2),xyzp(:,3),tri);
     xyzc1 = sum(repmat(at,1,3).*cgt)./sum(at);
%
     if iplt
       plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k.-', ...
             'LineWidth',1);

       plot3(xyzc1(:,1),xyzc1(:,2),xyzc1(:,3),'bs', ...
             'LineWidth',1);
     end
%
   end
%
% Posterior PCL
%
   if strcmp(roinaml,'postaxis')||strcmp(roinaml,'post_axis')
%
% Plot Posterior PCL
%
     if nslice>2
       warning([' *** WARNING in tibia_cs8:  Too many POSTAXIS ', ...
                'slices in ' faxial '!']);
     end
%
     nsz = zeros(2,1);
     xyz_post = zeros(2,3);
     for n = 1:2
         xyz = dat{n};
         nsz(n) = size(xyz,1);
         if iplt
           plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k^-', ...
                 'MarkerSize',8,'LineWidth',1);
         end
%
% Find Most Posterior Points on the Posterior Tibia Plateau at the PCL
%
         [~,id] = max(xyz(:,2));
         xyz_post(n,:) = xyz(id,:);
     end
%
% Combine Both Posterior Sides
%
     ns = sum(nsz);     % Total number of posterior points
     idp{2} = nsz(1)+1:ns;             % Index to second side
     idp{1} = 1:nsz(1);                % Index to first side
%
     xyz = cell2mat(dat);
     xyzo = mean(xyz);
     xyzro = xyz-repmat(xyzo,ns,1);    % Center data
%
% Loop until Y-Axis Changes are Less Than One (1) Degree
%
     while dang>ptol
%
          if leg==1
            if xyz_post(1,1)>xyz_post(2,1)
              xme(1) = xyz_post(1,1);
              yme(1) = xyz_post(1,2);
              zme(1) = xyz_post(1,3);
              xme(2) = xyz_post(2,1);
              yme(2) = xyz_post(2,2);
              zme(2) = xyz_post(2,3);
            else
              xme(1) = xyz_post(2,1);
              yme(1) = xyz_post(2,2);
              zme(1) = xyz_post(2,3);
              xme(2) = xyz_post(1,1);
              yme(2) = xyz_post(1,2);
              zme(2) = xyz_post(1,3);
            end
          else
            if xyz_post(1,1)<xyz_post(2,1)
              xme(1) = xyz_post(1,1);
              yme(1) = xyz_post(1,2);
              zme(1) = xyz_post(1,3);
              xme(2) = xyz_post(2,1);
              yme(2) = xyz_post(2,2);
              zme(2) = xyz_post(2,3);
            else
              xme(1) = xyz_post(2,1);
              yme(1) = xyz_post(2,2);
              zme(1) = xyz_post(2,3);
              xme(2) = xyz_post(1,1);
              yme(2) = xyz_post(1,2);
              zme(2) = xyz_post(1,3);
            end
          end
%
% Get Y-Axis
%
          yax = [xme(2)-xme(1) yme(2)-yme(1) zme(2)-zme(1)]';
%
% Get Rotation Angle and Rotate about the Z-Axis
%
          ang(1) = ang(2);
          ang(2) = -atan(yax(2)/yax(1));
          dang = diff(abs(ang));
%
          r = rotxyz([0 0 ang(2)]);
          xyzr = xyzro*r';
%
% Find Most Posterior Points on the Posterior Tibia Plateau at the PCL
%
          ids = zeros(2,1);
          for n = 1:2
             [~,id] = max(xyzr(idp{n},2));
             ids(n) = idp{n}(id);
             xyz_post(n,:) = xyz(ids(n),:);
          end
     end
%
% Plot Posterior Points and Line
%
     if iplt
       plot3(xyz_post(:,1),xyz_post(:,2),xyz_post(:,3), ...
             'bo','LineWidth',2);
     end
   end
%
% Get Proximal Tibia Section
%
   if strcmp(roinaml,'prox_tibia')
%
     if nslice>1
       warning([' *** WARNING in tibia_cs8:  Too many PROX_TIBIA ', ...
                'slices in ' faxial '!']);
     end
%
     xyz = dat{1};                     % Use only the first digitization
     npts = size(xyz,1);
%
     xyzm = mean(xyz);
     xyzp = [xyz; xyzm];
%
     if iplt
       plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k.-', ... 
             'LineWidth',1);
     end
%
% Get and Plot Centroid of Proximal Tibia Section
%
     tri = [repmat(npts+1,npts-1,1) (1:npts-1)' (2:npts)'];
     tri = [tri; [npts+1 npts 1]];
%
     [at,cgt] = tri_area(xyzp(:,1),xyzp(:,2),xyzp(:,3),tri);
     xyzc2 = sum(repmat(at,1,3).*cgt)./sum(at);
     maxheight=xyzc2(3);
     if iplt
       plot3(xyzc2(1),xyzc2(2),xyzc2(3),'bs', ...
             'LineWidth',1);
     end
   end
%
end
%
% Draw a Line Connecting Posterior Tibia Plateau Points
%
if iplt
  plot3(xme,yme,zme,'b-','LineWidth',2);
end
%
% Get Coordinate System
%
rvec = xyzc2 - xyzc1;   % Vertical Z axis
%
% Get Coordinate System Origin
%
xyzm = [xme yme zme];
xyzo = li_clos([xyzc1;xyzc2],xyzm);
xc = xyzo(1);
yc = xyzo(2);
zc = xyzo(3);
%
% Get Coordinate System
%
yax = [xme(2)-xme(1) yme(2)-yme(1) zme(2)-zme(1)]';
yax = yax./norm(yax);
zax = sign(rvec(3))*rvec';
zax = zax./norm(zax);
xax = cross(yax,zax);
xax = xax./norm(xax);
yax = cross(zax,xax);
yax = yax./norm(yax);   % Not necessary?
%
% Plot Coordinate System
%
if iplt
  h = quiver3(xc,yc,zc,xax(1),xax(2),xax(3),50,'k');
  set(h,'LineWidth',2);
  h = quiver3(xc,yc,zc,yax(1),yax(2),yax(3),50,'k');
  set(h,'LineWidth',2);
  h = quiver3(xc,yc,zc,zax(1),zax(2),zax(3),50,'k');
  set(h,'LineWidth',2);
  axis equal;           % Make axes equal scale
end
%
% Return Axis Vectors in Matrices
%
xyzax = [xax yax zax];
%
% Calculate Aspect Ratio Based on Maximum Width of the Proximal Tibia
% in the Tibial Coordinate System and Height Between the Distal and
% Proximal Tibia Outlines
%
npts = size(xyzp,1);
xyzp = xyzp - repmat(xyzo,npts,1);
xyzp = xyzp*xyzax;
widt = abs(max(xyzp(:,1))-min(xyzp(:,1)));
minheight = xyzc1(3);
maxheight = xyzc2(3);
height = maxheight-minheight;
aspect = widt/height;
%
if nargout==0
  fprintf(1,'\n Aspect ratio = %5.3f\n\n');
end
%
return