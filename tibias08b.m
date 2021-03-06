%#######################################################################
%
%                      * TIBIAS 08 Bone Program *
%
%          Program to read and plot digitized points from a tibia in the
%     tibia coordinate system.  Plots medial and lateral compartment
%     subchondral bone data for user selected knee.  The tibial
%     coordinate system for the knee is saved to a Matlab MAT file.
%
%     NOTES:  1.  The names of the regions of interest (ROIs) must
%             follow a specific convention.
%
%             2.  See tibia_cs8.m for more information on the tibia
%             coordinate system.
%
%             3. The Matlab M-files fix_pts.m, li_clos.m, mk_tri6.m,
%             plane_fit.m, rd_roi4.m, rotxyz.m, tibia_cs8.m,
%             tri_area.m, tri_fix2.m and tri_norm.m must be in the
%             current directory or path.
%
%             4.  This file outputs a PostScript file, tibias08b.ps,
%             with plots of the tibia coordinate system and the
%             subchondral bone into the directory of the knee data.
%
%             5.  The subchondral bone data in the tibia coordinate
%             system and the coordinate transformation matrix and
%             origin are saved in the directory of the knee data in a
%             Matlab MAT file, kid_tibiaCS.mat, where kid is the knee
%             identifier (xxx_R/L, xxx is the knee number and either R
%             or L for the right or left knee).
%
%             6.  The aspect ratio of the tibia width to the tibia
%             height is out to a CSV file, kid_AspectRatio.csv.
%
%             7.  This is an updated version of tibias07b.m.
%
%     07-May-2019 * Mack Gardner-Morse
%

%#######################################################################
%
% Control Variables for Finding Duplicate Points
%
tol = 0.2;              % Minimum distance between distinct points
iflag = true;           % Print message if duplicates found
%
% MRI Slice Thickness
%
sl_thick = 0.3125;      % MRI slice thickness in mm
%
% Output PS and MAT File Names
%
psnam = '_tibias08b.ps';
mnam = '_tibiaCS.mat';
%
% Compartment Colors and Labels
%
tclrs = [0 0 0.7; 0 0.5 0];            % Deep blue and dark green
%
tcmpt = ['lateral'
         'medial '];
%
% Get Sagittal Bone CSV File Name
%
if ~exist('fnam','var')
  [fnam,pnam] = uigetfile({'*_SAG_TIB*.csv', ...
                  'Sagittal tibial bone CSV files'}, ...
                  'Please Select Tibial Sagittal Bone CSV Files');
end
%
if isequal(fnam,0)
  return;
end
%
% Get Knee ID
%
kid = fnam(1:5);
%
% Check for Axial Plane File
%
fnama = [kid '_AX_*.csv'];
fnamas = dir(fullfile(pnam,fnama));
%
if isempty(fnamas)
  error(' *** ERROR in tibias08b:  Unable to find axial file!');
end
%
fnamas = {fnamas.name}';     % Get axial file names in a cell array
%
% Find Tibia Axial File
%
itib = ~contains(upper(fnamas),'FEM'); % Not a femur axial file
fnama = char(fnamas(itib,:));
%
if size(fnama,1)~=1
  error([' *** ERROR in tibias08b:  Unable to find an unique ', ....
        'axial file!']);
end
%
% Get Leg
%
ileg = strcmpi(kid(5),'R');
%
% Setup Tibial Coordinate System Figure
%
hf1 = figure;
orient tall;
view(3);
hold on;
%
% Get Tibial Coordinate System, Proximal Tibia Outline (PTO) and Aspect
% Ratio from the PTO
%
% NOTE:  PTO is in the tibia coordinate system.
%
[xyzc,xyzr,aspect,widt,height,xyzpto] = tibia_cs8(fullfile(pnam,fnama), ...
                                                  ileg,true);
%
% Finish Plot
%
title({kid; 'Tibia Coordinate System'},'FontSize',16, ...
      'FontWeight','bold', 'Interpreter','none');
axis equal;
%
psnam = [kid psnam];
psnam = fullfile(pnam,psnam);
%
print('-dpsc2','-r300','-bestfit',psnam);
%
% Output Aspect Ratio to Screen and a CSV File
%
diff_ht = widt-height;  % Width and height should be equal
nchange_sl = round(diff/sl_thick);     % Number of slices to move
%
fprintf(1,'\n\n Aspect Ratio = %.3f\n',aspect)
if nchange_sl>0
  fprintf(1,'   Go Distally %.0f Slices\n\n',abs(nchange_sl));
elseif nchange_sl<0
  fprintf(1,'   Go Proximally %.0f Slices\n\n',abs(nchange_sl));
else
  fprintf(1,'   Aspect Ratio is Close!\n\n');
end
%
fid = fopen(fullfile(pnam,[kid '_AspectRatio.csv']),'w');
%
fprintf(fid,'"kid","AR","Width","Height"\n');
fprintf(fid,',,"(mm)","(mm)"\n');
fprintf(fid,'%s,%g,%g,%g\n',fnama,aspect,widt,height);
%
fclose(fid);
%
% Get Raw Tibia Bone Data
%
roi = rd_roi4(fullfile(pnam,fnam));
%
roinams = upper(char(roi.name));
idc(2) = find(strcmp('MTCS',cellstr(roinams)));  % Medial compartment
idc(1) = find(strcmp('LTCS',cellstr(roinams)));  % Lateral compartment
%
% Raw Data Figure
%
hf2 = figure;
orient landscape;
view(3);
hold on;
xlabel('X (mm)','FontSize',12,'FontWeight','bold');
ylabel('Y (mm)','FontSize',12,'FontWeight','bold');
zlabel('Z (mm)','FontSize',12,'FontWeight','bold');
title({[kid ' - MRI CS']; 'Blue - Lateral, Green - Medial'}, ...
      'FontSize',16,'FontWeight','bold','Interpreter','none');
%
% Transformed Data Figure
%
hf3 = figure;
orient landscape;
view(3);
hold on;
xlabel('AP (mm)','FontSize',12,'FontWeight','bold');
ylabel('Lateral (mm)','FontSize',12,'FontWeight','bold');
zlabel('Superior (mm)','FontSize',12,'FontWeight','bold');
title([kid ' - Tibial CS'], ...
      'FontSize',16,'FontWeight','bold','Interpreter','none');
%
% Loop through Compartments
%
for l = 1:2
%
% Get Data
%
   dat1 = roi(idc(l)).data';
   nsl1 = size(dat1,1);              % Number of slices
%
% Check for Duplicates, Direction of Digitization and Transform Data to
% Tibia Coordinate System
%
   datt = cell(nsl1,1); % Transformed data
%
   for n = 1:nsl1
      xyz = dat1{n};
      if isempty(xyz)
        cmprt = deblank(tcmpt(l,:));
        error([' *** ERROR in tibias08b:  No coordinates for ', ...
               'slice ' int2str(n) ' in ' cmprt ' compartment!']);
      end
      xyz = fix_pts(xyz,tol,iflag);    % Check for duplicates
      npts = size(xyz,1);
      [~,imx] = max(xyz(:,2));
      [~,imn] = min(xyz(:,2));
      if imn<imx
        xyz = xyz(1:npts,:);           % Anterior to posterior
      else
         xyz = xyz(npts:-1:1,:);       % Anterior to posterior
      end
      xyzt = xyz-repmat(xyzc,npts,1); % Center data
      xyzt = xyzt*xyzr;               % Transform data to tibia CS
      datt{n} = xyzt;
% 
% Plot Raw Data
%
      figure(hf2);
%
      plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k.-','Color',tclrs(l,:), ...
            'MarkerSize',8,'LineWidth',1);
      if xyz(1,2)<0
        text(xyz(1,1),xyz(1,2)-0.5,xyz(1,3),int2str(n), ...
             'Color','k','FontSize',10);
      else
        text(xyz(1,1),xyz(1,2)+0.5,xyz(1,3),int2str(n), ...
             'Color','k','FontSize',10);
      end
% 
% Plot Transformed Data
%
      figure(hf3);
%
      plot3(xyzt(:,1),xyzt(:,2),xyzt(:,3),'k.-','Color',tclrs(l,:), ...
            'MarkerSize',8,'LineWidth',1);
      if xyzt(1,1)<0
        text(xyzt(1,1)-0.5,xyzt(1,2),xyzt(1,3),int2str(n), ...
             'Color','k','FontSize',10);
      else
        text(xyzt(1,1)+0.5,xyzt(1,2),xyzt(1,3),int2str(n), ...
             'Color','k','FontSize',10);
      end
   end
%
% Get Surface Triangulation
%
   trit = mk_tri6(datt);
   xyzt = cell2mat(datt);
   trit = tri_fix2(trit,xyzt); %#ok<NASGU>
%
% Save Data into Compartment Specific Variables
%
   eval(['dat' tcmpt(l,1) ' = datt;']);
   eval(['tri' tcmpt(l,1) ' = trit;']);
   eval(['xyz' tcmpt(l,1) ' = xyzt;']);
%
   clear datt trit xyzt;
%
end
%
% Save Plots
%
print(hf2,'-dpsc2','-r300','-bestfit','-append',psnam);
print(hf3,'-dpsc2','-r300','-bestfit','-append',psnam);
%
% Plot Transformed Bone Surface Data
% 
hf4 = figure;
orient landscape;
%
trisurf(tril,xyzl(:,1),xyzl(:,2),xyzl(:,3),'FaceColor','interp', ...
        'EdgeColor',tclrs(1,:),'LineWidth',1);
hold on;
trisurf(trim,xyzm(:,1),xyzm(:,2),xyzm(:,3),'FaceColor','interp', ...
        'EdgeColor',tclrs(2,:),'LineWidth',1);
axis equal;
xlabel('AP (mm)','FontSize',12,'FontWeight','bold');
ylabel('Lateral (mm)','FontSize',12,'FontWeight','bold');
zlabel('Superior (mm)','FontSize',12,'FontWeight','bold');
title([kid ' - Tibial CS'], ...
      'FontSize',16,'FontWeight','bold','Interpreter','none');
%
print('-dpsc2','-r300','-opengl','-bestfit','-append',psnam);
%
% Save Data into a Matlab MAT File for Further Processing
%
mnam = [kid mnam];
mnam = fullfile(pnam,mnam);
%
save(mnam,'datl','datm','kid','ileg','tril','trim','xyzc','xyzr', ...
          'xyzl','xyzm','xyzpto');
%
return