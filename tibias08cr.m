%#######################################################################
%
%              * TIBIAS 08 Cartilage Reliability Program *
%
%          Program to read, transform to the bony tibia coordinate
%     system and plot digitized cartilage points from a tibia.  Plots
%     lateral compartment cartilage data for user selected knee.  Plots
%     both coronal and sagittal data.  Modified version of tibias08c.m
%     for reliability study.
%
%     NOTES:  1.  Reliability data must be in a subdirectory of the
%             knee data folder.
%
%             2.  The names of the regions of interest (ROIs) must
%             follow a specific convention.
%
%             3. The Matlab M-files fix_pts.m, li_clos.m, mk_tri6.m,
%             plane_fit.m, rd_roi4.m, rotxyz.m, tri_area.m, tri_fix2.m
%             and tri_norm.m must be in the current directory or path.
%
%             4.  The tibia coordinate system with the coordinate
%             transformation matrix and origin are read from the
%             directory of the knee data from the Matlab MAT file,
%             kid_tibiaCS.mat, where kid is the knee identifier
%             (xxx_R/L, xxx is the knee number and either R or L
%             for the right or left knee).
%
%             5.  This M-file outputs a PostScript file,
%             kid_tibias08cr.ps, with plots of the tibial cartilage
%             into the reliability subdirectory of the knee data.
%
%             6.  The transformed data and mesh are saved in a Matlab
%             MAT file kid_tibiaCartR.mat into the reliability
%             subdirectory.
%
%     05-Dec-2019 * Mack Gardner-Morse
%

%#######################################################################
%
% Control Variables for Finding Duplicate Points
%
tol = 0.2;              % Minimum distance between distinct points
iflag = true;           % Print message if duplicates found
%
% Tibial Coordinate and Output PS and MAT File Names
%
csnam = '_tibiaCS.mat';                % Tibia coordinate system file
psnam = '_tibias08cr.ps';              % Output PS file
mnam = '_tibiaCartR.mat';              % Output MAT file
%
% Digitization Planes (Coronal and Sagittal)
%
idig = ['c'; 's'];      % Coronal and sagittal
%
% Compartment Colors and Labels
%
tclrs = [0 0 0.7; 0 0.5 0];            % Deep blue and dark green
%
tcmpt = ['lateral'
         'medial '];
%
% Get Sagittal Cartilage CSV File Name
%
if ~exist('fnams','var')
  [fnams,pnam] = uigetfile({'*_SAGAR_TIB*.csv', ...
                   'Sagittal tibial cartilage CSV files'}, ...
                   'Please Select Tibial Sagittal Cartilage CSV Files');
end
%
if isequal(fnams,0)
  return;
end
%
% Get Knee ID
%
kid = fnams(1:5);
%
% Check for Coronal Plane File
%
fnamc = [kid '_CORAR_TIB*.csv'];
fnamc = dir(fullfile(pnam,fnamc));
%
if isempty(fnamc)
  error(' *** ERROR in tibias08cr:  Unable to find coronal file!');
end
%
if size(fnamc,1)~=1
  error([' *** ERROR in tibias08cr:  Unable to find an unique ', ....
        'coronal file!']);
end
%
fnamc = fnamc.name;
%
% Get Leg
%
ileg = strcmpi(kid(5),'R');
%
% Get Tibial Coordinate System
%
idp = find(pnam(1:end-1)==filesep);    % Skip any file separator at end of path
idp = idp(end);         % Remove subdirectory from path
pnamb = pnam(1:idp);    % Path to bone data
%
tcs = load(fullfile(pnamb,[kid csnam]),'xyzc','xyzr');
xyzc = tcs.xyzc;        % Origin
xyzr = tcs.xyzr;        % Rotation matrix
%
% Full PS File Name
%
psnam = [kid psnam];
psnam = fullfile(pnam,psnam);
%
% Loop through Digitizations (Coronal [k==1] and Sagittal [k==2])
%
for k = 1:2             % Coronal = 1 and sagittal = 2
%
% Get Raw Tibia Cartilage Data
%
   if k==1
     roi = rd_roi4(fullfile(pnam,fnamc));
     roinams = upper(char(roi.name));
     idc(2) = find(strcmp('MACC',cellstr(roinams)));  % Medial compartment
     idc(1) = find(strcmp('LACC',cellstr(roinams)));  % Lateral compartment
     dtxt = 'Coronal Digitization';
   else
     roi = rd_roi4(fullfile(pnam,fnams));
     roinams = upper(char(roi.name));
     idc(2) = find(strcmp('MACS',cellstr(roinams)));  % Medial compartment
     idc(1) = find(strcmp('LACS',cellstr(roinams)));  % Lateral compartment
     dtxt = 'Sagittal Digitization';
   end
%
% Raw Data Figure
%
   hf1 = figure;
   orient landscape;
   view(3);
   hold on;
   xlabel('X (mm)','FontSize',12,'FontWeight','bold');
   ylabel('Y (mm)','FontSize',12,'FontWeight','bold');
   zlabel('Z (mm)','FontSize',12,'FontWeight','bold');
   title({[kid ' - MRI CS']; dtxt; 'Lateral Only'}, ...
         'FontSize',16,'FontWeight','bold', ...
         'Interpreter','none');
%
% Transformed Data Figure
%
   hf2 = figure;
   orient landscape;
   view(3);
   hold on;
   xlabel('AP (mm)','FontSize',12,'FontWeight','bold');
   ylabel('Lateral (mm)','FontSize',12,'FontWeight','bold');
   zlabel('Superior (mm)','FontSize',12,'FontWeight','bold');
   title({[kid ' - Tibial CS']; dtxt}, ...
         'FontSize',16,'FontWeight','bold','Interpreter','none');
%
% Loop through Compartments (Lateral [l==1] and Medial [l==2])
%
   for l = 1:1          % Lateral only
%    for l = 1:2          % Lateral = 1 and medial = 2
%
% Get Data
%
      dat1 = roi(idc(l)).data';
      nsl1 = size(dat1,1);             % Number of slices
%
% Check for Duplicates, Direction of Digitization and Transform Data to
% Tibia Coordinate System
%
      datt = cell(nsl1,1);             % Transformed data
%
      for n = 1:nsl1
         xyz = dat1{n};
         if isempty(xyz)
           cmprt = deblank(tcmpt(l,:));
           error([' *** ERROR in tibias08c:  No coordinates for ', ...
                  dtxt ' slice ' int2str(n) ' in ' cmprt, ...
                  'compartment!']);
         end
         xyz = fix_pts(xyz,tol,iflag); % Check for duplicates
         npts = size(xyz,1);
         [~,imx] = max(xyz(:,2));
         [~,imn] = min(xyz(:,2));
         if imn<imx
           xyz = xyz(1:npts,:);        % Anterior to posterior
         else
            xyz = xyz(npts:-1:1,:);    % Anterior to posterior
         end
         xyzt = xyz-repmat(xyzc,npts,1);    % Center data
         xyzt = xyzt*xyzr;                  % Transform data to tibia CS
         datt{n} = xyzt;
% 
% Plot Raw Data
%
         figure(hf1);
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
         figure(hf2);
%
         plot3(xyzt(:,1),xyzt(:,2),xyzt(:,3),'k.-','Color', ...
               tclrs(l,:),'MarkerSize',8,'LineWidth',1);
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
      trit = tri_fix2(trit,xyzt); 
%
% Save Data into Compartment Specific Variables
%
      eval(['dat' tcmpt(l,1) idig(k) ' = datt;']);
      eval(['tri' tcmpt(l,1) idig(k) ' = trit;']);
      eval(['xyz' tcmpt(l,1) idig(k) ' = xyzt;']);
%
      clear datt trit xyzt;
   end
%
% Save Plots
%
   print(hf1,'-dpsc2','-r300','-bestfit','-append',psnam);
   print(hf2,'-dpsc2','-r300','-bestfit','-append',psnam);
%
% Plot Transformed Cartilage Surface Data
% 
   hf3 = figure;
   orient landscape;
%
   if k==1
%
     trisurf(trilc,xyzlc(:,1),xyzlc(:,2),xyzlc(:,3),'FaceColor', ...
             'interp','EdgeColor',tclrs(1,:),'LineWidth',1);
     hold on;
%      trisurf(trimc,xyzmc(:,1),xyzmc(:,2),xyzmc(:,3),'FaceColor', ...
%              'interp','EdgeColor',tclrs(2,:),'LineWidth',1);
   else
     trisurf(trils,xyzls(:,1),xyzls(:,2),xyzls(:,3),'FaceColor', ...
             'interp','EdgeColor',tclrs(1,:),'LineWidth',1);
     hold on;
%      trisurf(trims,xyzms(:,1),xyzms(:,2),xyzms(:,3),'FaceColor', ...
%              'interp','EdgeColor',tclrs(2,:),'LineWidth',1);
   end
   axis equal;
   xlabel('AP (mm)','FontSize',12,'FontWeight','bold');
   ylabel('Lateral (mm)','FontSize',12,'FontWeight','bold');
   zlabel('Superior (mm)','FontSize',12,'FontWeight','bold');
   title({[kid ' - Tibial CS']; dtxt}, ...
         'FontSize',16,'FontWeight','bold','Interpreter','none');
%
   print('-dpsc2','-r300','-opengl','-bestfit','-append',psnam);
%
end
%
% Save Data into a Matlab MAT File for Further Processing
%
mnam = [kid mnam];
mnam = fullfile(pnam,mnam);
%
save(mnam,'datlc','datls','kid','ileg','trilc','trils','xyzlc','xyzls');
%
return