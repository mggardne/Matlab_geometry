%#######################################################################
%
%          * Tibia CARilage Thickness 8 Reliability Program *
%
%          M-File which runs through the left or right tibia
%     within selected subject directories and defines a 1 mm by 1 mm
%     grid, scales and projects the grid onto the bone surface mesh,
%     calculates the intersection of the bone normals with the
%     cartilage surface mesh and calculates cartilage thicknesses.  The
%     same grid is used for both the sagittal and coronal cartilage
%     digitizations, so similar points can be averaged across the
%     digitizations.  Modified version of tcart08.m for reliability
%     study.
%
%     NOTES:  1.  The M-files car_thk6.m, gridproj.m, line_fit.m,
%             nod2tri.m, nod_norm.m, pts2lin.m, tri_fix2.m,
%             tri_norm.m, tsect4.m and xprod.m must be in the current
%             path or directory.
%
%             2.  Bone and cartilage MAT files must already exist in
%             the selected subject directories.
%
%             3.  This program produces a MAT file with the scaling for
%             the tibias within the selected directories: tgrid08_*.mat,
%             where * is the number of selected subject directories.
%             Final thicknesses are saved in the MAT file:
%             tcart08_*.mat.  The MAT files are saved in the directory:
%             TibiaCartThkR_*_ddmmmyyyy, where * is the number of
%             selected subject directories and ddmmmyyyy is the date in
%             day, month (first three letters) and year format.
%
%             4.  This M-file outputs PostScript files,
%             tcart08_coronal.ps, tcart08_sagittal.ps, 
%             tcart08_coronal_thk.ps and tcart08_sagittal_thk.ps with
%             plots of the tibial bony and cartilage surfaces and tibial
%             cartilage thicknesses in the results directory.
%
%     05-Dec-2019 * Mack Gardner-Morse
%

%#######################################################################
%
% Plot Parameter
%
% iplt = false;           % No plots
iplt = true;            % Plots
%
% iprt = false;           % No printing of plots
iprt = true;            % Print plots
%
% Get All Subject Subdirectories
%
ddir = fullfile('..','CACL MRI Master');
sd = dir(ddir);         % All subdirectories
snams = {sd.name}';
id = startsWith(snams,{'.'; '..'});    % Current and parent directories
snams = snams(~id);     % Just subject directories
%
% Select Subject Subdirectories for Analysis
%
idx = listdlg('ListString',snams,'ListSize',[200 600],'Name', ...
              'Tibias','PromptString', ...
              {'    Select all subdirectories for'; ...
              'calculating cartilage thicknesses.'});
%
if isempty(idx)
  return;
end
%
snams = char(snams(idx));    % Make subject names into character array
%
ns = size(snams,1);     % Number of subjects
nss = int2str(ns);      % # of subjects for results directory and file
%
% Get Results Directory Name
%
dstr = date;            % Get date as a string
dstr = strrep(dstr,'-','');  % Remove dashes
%
rdir = fullfile(ddir,['TibiaCartThkR_' nss '_' dstr]);
%
% Initialize Coordinate and Bone Arrays
%
ilegs = false(ns,1);    % 0 (false) for left/1 (true) for right
kids = repmat(blanks(ns)',1,5);        % Knee IDs
tribl = cell(ns,1);     % Sagittal bone lateral triangular mesh connectivity
xyzbl = cell(ns,1);     % Sagittal lateral bone point coordinates
xyzmnl = zeros(ns,3);   % Minimum individual lateral tibia coordinates
xyzmxl = zeros(ns,3);   % Maximum individual lateral tibia coordinates
xyzs = zeros(ns,3);     % Range of proximal tibia outline
%
% Loop through Subjects (Tibias)
%
for k = 1:ns
%
% Get Subject Number
%
   snam = snams(k,:);   % Subject directory
%
   bfnam = fullfile(ddir,snam,[snam(1:3) '_*_tibiaCS.mat']);    % Bone MAT file name
   d = dir(bfnam);
   if isempty(d)
     error([' *** ERROR in tcart08r:  Not able to find tibia bone ', ...
            'file:  ',bfnam,'\n']);
   else
     bfnam = fullfile(ddir,snam,d.name);    % Bone MAT file name
   end
%
% Get Bone Data
%
   bone = load(bfnam);
   if isfield(bone,'xyzpto')
     ilegs(k) = bone.ileg;
     kids(k,:) = bone.kid;
     tribl{k} = bone.tril;
     xyzbl{k} = bone.xyzl;
     xyzpt = bone.xyzpto;
   else
     error([' *** ERROR in tcart08r:  Not able to find proximal ', ...
            'tibial outline.  Please run tibias08b again.']);
   end
%
% Reverse X-Axis for Right Knees
%
   if ilegs(k)
     xyzbl{k}(:,1) = -xyzbl{k}(:,1);
   end
%
% Get Minimum and Maximum Coordinates
%
   xyzmnl(k,:) = min(xyzbl{k});
   xyzmxl(k,:) = max(xyzbl{k});
%
% Get Range of Proximal Tibia Outline
%
   xyzs(k,:) = max(xyzpt)-min(xyzpt);
%
end
%
% Generate Scaling (Based on Proximal Tibia Outline) and Uniform Grid
%
xyzs_avg = mean(xyzs);
sc = xyzs./repmat(xyzs_avg,ns,1);
scx = sc(:,1);          % X scale
scy = sc(:,2);          % Y scale
%
% Get Range and Make a Uniform Rectangular Grid in X and Y for the
% Lateral Compartment
%
x = (floor(xyzmnl(1)):ceil(xyzmxl(1)))';    % X range in 1 mm increments
nxl = size(x,1);
y = (floor(xyzmnl(2)):ceil(xyzmxl(2)))';    % Y range in 1 mm increments
nyl = size(y,1);
[X,Y] = meshgrid(x,y);  % Rectangle grid of square quadrilaterals
xql = X(:);             % Grid points as a vector
yql = Y(:);             % Grid points as a vector
nl = size(xql,1);       % Number of grid points in lateral compartment
%
quadl = (1:nyl-1)';
quadl = [quadl quadl+nyl quadl+nyl+1 quadl+1];
quadl = repmat(quadl,nxl-1,1);
addcol = repmat(nyl*(0:nxl-2),(nyl-1)*4,1);
addcol = reshape(addcol,4,(nxl-1)*(nyl-1))';
quadl = quadl+addcol;
%
nel = size(quadl,1);    % Number of lateral quadrilateral elements
%
% Triangular Grids
%
trigl = [quadl(:,1:3) quadl(:,1) quadl(:,3:4)]';
trigl = reshape(trigl,3,2*nel)';
%
% Save Analysis Grids
%
gnam = fullfile(rdir,['tgrid08_' nss '.mat']);   % Grid MAT file
%
if exist(rdir,'dir')
  ierr = menu('Overwrite Previous Results','Yes','No')-1;
  if ierr
    return;
  end
else
  mkdir(rdir);
end
%
save(gnam,'nel','nl','nxl','nyl','ns','quadl','sc','scx','scy', ...
     'trigl','xql','yql','xyzmnl','xyzmxl');
%
% Initialize Cartilage Arrays
%
trilcs = cell(ns,1);    % Lateral coronal triangular mesh
trilss = cell(ns,1);    % Lateral sagittal triangular mesh
xyzlcs = cell(ns,1);    % Lateral coronal coordinates
xyzlss = cell(ns,1);    % Lateral sagittal coordinates
%
% Initialize Cartilage Thickness Arrays
%
cthklc = zeros(nl,ns);       % Lateral coronal cartilage thicknesses
xyzilc = zeros(nl,3,ns);     % Lateral coronal cartilage intersection
%
cthkls = zeros(nl,ns);       % Lateral sagittal cartilage thicknesses
xyzils = zeros(nl,3,ns);     % Lateral sagittal cartilage intersection
%
zgl = zeros(nl,ns);     % Most superior Z coordinates on bony grid surface
%
% Loop through Subjects (Tibias)
%
for k = 1:ns
%
% Get Subject Number
%
   snam = snams(k,:);   % Subject directory
   kid = kids(k,:);     % Knee ID
   pnam = fullfile(ddir,snam);
%
% Get Reliability Subdirectory
%
   ddirr = dir(pnam);
   ddirr = {ddirr([ddirr.isdir]).name}';
   idr = startsWith(ddirr,'rel','IgnoreCase',true);
   if any(idr)
     ddirr = ddirr{idr};
     if size(ddirr,1)==1
       pnamr = fullfile(pnam,ddirr);
     else
       warning([' *** WARNING in tcart08r:  Unable to find ', ...
                'reliability subdirectory for ',snam,'!']);
       break;
     end
   else
     warning([' *** WARNING in tcart08r:  Unable to find ', ...
              'reliability subdirectory for ',snam,'!']);
     break;
   end
%
% Get Cartilage Data
%
   cfnam = fullfile(pnamr,[kid '_tibiaCartR.mat']);   % Cartilage MAT file name
   if ~exist(cfnam,'file')
     error([' *** ERROR in tcart08r:  Not able to find tibia ', ...
            'cartilage file:  ',cfnam,'\n']);
   end
   cart = load(cfnam);
%
% Get and Save Cartilage Data
%
   trilcs{k} = cart.trilc;
   trilss{k} = cart.trils;
   xyzlcs{k} = cart.xyzlc;
   xyzlss{k} = cart.xyzls;
%
% Reverse X-Axis for Right Knees
%
   if ilegs(k)
     xyzlcs{k}(:,1) = -xyzlcs{k}(:,1);
     xyzlss{k}(:,1) = -xyzlss{k}(:,1);
   end
%
% Scale Grid
%
   xgl = xql*scx(k);
   ygl = yql*scy(k);
%    xgl = xql;           % No scaling
%    ygl = yql;           % No scaling
%
% Get Bony Tibia Data
%
   xyzl = xyzbl{k};     % Lateral coordinates
%
   tril = tribl{k};     % Lateral triangular mesh connectivity
%
% Project Grid onto Bony Surface Mesh
%
   zgl(:,k) = gridproj(tril,xyzl,xgl,ygl,1);
%
   xyzgl = [xgl ygl zgl(:,k)];
%
% Improve the Mesh
%
   trigli = tri_fix2(trigl,xyzgl);
%
% Plot Cartilage and Bone Surfaces
%
   if iplt
%
% Plot Coronal Cartilage and Gridded Bone Surfaces
%
     if exist('hf1','var')
       figure(hf1);
       clf;
     else
       hf1 = figure('Position',[1 53 1920 952]);
     end
     hb1l = trimesh(trigli,xgl,ygl,zgl(:,k), ...
                    'LineWidth',0.25,'FaceColor','none', ...
                    'EdgeColor','k');
     hold on;
%
     hc1l = trimesh(trilcs{k},xyzlcs{k}(:,1),xyzlcs{k}(:,2), ...
                    xyzlcs{k}(:,3),'LineWidth',0.25, ...
                    'FaceColor','none','EdgeColor','b');
%
     xlabel('X','FontSize',12,'FontWeight','bold');
     ylabel('Y','FontSize',12,'FontWeight','bold');
     zlabel('Z','FontSize',12,'FontWeight','bold');
     title({kid; 'Coronal Cartilage and Bone Meshes'}, ...
           'Interpreter','none','FontSize',16,'FontWeight','bold');
     view(-60,12);
     axis equal;
     if iprt
       psnam = fullfile(rdir,'tcart08_coronal.ps');
       orient landscape;
       if k==1
         print('-dpsc2','-r300','-fillpage',psnam);
       else
         print('-dpsc2','-r300','-fillpage','-append',psnam);
       end
     end
%
% Plot Sagittal Cartilage and Gridded Bone Surfaces
%
     if exist('hf2','var')
       figure(hf2);
       clf;
     else
       hf2 = figure('Position',[1 53 1920 952]);
     end
     hb2l = trimesh(trigli,xgl,ygl,zgl(:,k), ...
                    'LineWidth',0.25,'FaceColor','none', ...
                    'EdgeColor','k');
     hold on;
%
     hc2l = trimesh(trilss{k},xyzlss{k}(:,1),xyzlss{k}(:,2), ...
                    xyzlss{k}(:,3),'LineWidth',0.25, ...
                    'FaceColor','none','EdgeColor','b');
     xlabel('X','FontSize',12,'FontWeight','bold');
     ylabel('Y','FontSize',12,'FontWeight','bold');
     zlabel('Z','FontSize',12,'FontWeight','bold');
     title({kid; 'Sagittal Cartilage and Bone Meshes'}, ...
           'Interpreter','none','FontSize',16,'FontWeight','bold');
     view(-60,12);
     axis equal;
     if iprt
       psnam = fullfile(rdir,'tcart08_sagittal.ps');
       orient landscape;
       if k==1
         print('-dpsc2','-r300','-fillpage',psnam);
       else
         print('-dpsc2','-r300','-fillpage','-append',psnam);
       end
     end
%      pause
   end
%
% Calculate and Save Coronal Cartilage Thicknesses
%
   [ct,bi] = car_thk6(trilcs{k},xyzlcs{k},trigli,xyzgl);   % Lateral
   cthklc(:,k) = ct;    % Save thicknesses
   xyzilc(:,:,k) = bi;  % Save cartilage intersection points
%
% Calculate and Save Sagittal Cartilage Thicknesses
%
   [ct,bi] = car_thk6(trilss{k},xyzlss{k},trigli,xyzgl);   % Lateral
   cthkls(:,k) = ct;    % Save thicknesses
   xyzils(:,:,k) = bi;  % Save cartilage intersection points
%
% Plot the Cartilage Thicknesses
%
   if iplt
%
% Plot Coronal Cartilage Thicknesses
%
     if exist('hf3','var')
       figure(hf3);
       clf;
     else
       hf3 = figure('Position',[1 53 1920 952]);
     end
     hb3l = trimesh(trigli,xgl,ygl,zgl(:,k), ...
                    'LineWidth',0.25,'FaceColor','none', ...
                    'EdgeColor','b');
     hold on;
%
     ctl = squeeze(cthklc(:,k));
     idx = find(~isnan(ctl));
     it = nod2tri(idx,trigli,2);
     hs3l = trisurf(trigli(it,:),xgl,ygl,zgl(:,k),ctl, ...
                    'FaceColor','interp', ...
                    'EdgeColor','b','LineWidth',0.25);
     colormap jet;
     view(-60,12);
     axis equal;
%
     ct = ctl;
     if min(ct)<0
       caxis([min(ct) max(ct)]);
     else
       caxis([0 max(ct)]);
     end
     hcb3 = colorbar;
     xlabel('X','FontSize',12,'FontWeight','bold');
     ylabel('Y','FontSize',12,'FontWeight','bold');
     zlabel('Z','FontSize',12,'FontWeight','bold');
     title({kid; 'Coronal Cartilage Thicknesses'},'Interpreter', ...
           'none','FontSize',16,'FontWeight','bold');
     if iprt
       psnam = fullfile(rdir,'tcart08_coronal_thk.ps');
       orient landscape;
       if k==1
         print('-dpsc2','-r300','-fillpage',psnam);
       else
         print('-dpsc2','-r300','-fillpage','-append',psnam);
       end
     end
%
% Plot Sagittal Cartilage Thicknesses
%
     if exist('hf4','var')
       figure(hf4);
       clf;
     else
       hf4 = figure('Position',[1 53 1920 952]);
     end
     hb4l = trimesh(trigli,xgl,ygl,zgl(:,k), ...
                    'LineWidth',0.25,'FaceColor','none', ...
                    'EdgeColor','b');
     hold on;
%
     ctl = squeeze(cthkls(:,k));
     idx = find(~isnan(ctl));
     it = nod2tri(idx,trigli,2);
     hs4l = trisurf(trigli(it,:),xgl,ygl,zgl(:,k),ctl, ...
                    'FaceColor','interp', ...
                    'EdgeColor','b','LineWidth',0.25);
     colormap jet;
     view(-60,12);
     axis equal;
%
     ct = ctl;
     if min(ct)<0
       caxis([min(ct) max(ct)]);
     else
       caxis([0 max(ct)]);
     end
     hcb4 = colorbar;
     xlabel('X','FontSize',12,'FontWeight','bold');
     ylabel('Y','FontSize',12,'FontWeight','bold');
     zlabel('Z','FontSize',12,'FontWeight','bold');
     title({kid; 'Sagittal Cartilage Thicknesses'}, ...
    'Interpreter','none','FontSize',16,'FontWeight','bold');
     if iprt
       psnam = fullfile(rdir,'tcart08_sagittal_thk.ps');
       orient landscape;
       if k==1
         print('-dpsc2','-r300','-fillpage',psnam);
       else
         print('-dpsc2','-r300','-fillpage','-append',psnam);
       end
     end
   end                  % End if iplt
%
end
%
% Save Tibia Data
%
rnam = fullfile(rdir,['tcart08_' nss '.mat']);    % Results MAT file
save(rnam,'cthklc','cthkls','ilegs','kids','snams', ...
     'tribl','trilcs','trilss','xyzbl','xyzilc','xyzils', ...
     'xyzlcs','xyzlss','zgl');
%
close all;              % Close cartilage thickness plots
%
return