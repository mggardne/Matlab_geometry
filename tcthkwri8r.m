%#######################################################################
%
%       * Tibia Cartilage THicKness WRIte 8 Reliability Program *
%
%          M-File which reads the tibia lateral compartment cartilage
%     thicknesses for all 21 subjects in the reliability study and
%     writes out the thicknesses to three MS-Excel spreadsheets:
%     tcthk_ant_lat.xlsx, tcthk_ctr_lat.xlsx and tcthk_pos_lat.xlsx.
%
%     NOTES:  1.  MAT files tgrid08_21.mat and tcart08_21.mat must be
%             in the two directories:
%             ..\CACL MRI Master\TibiaCartThk_21_13Dec2019
%             ..\CACL MRI Master\TibiaCartThkR_21_13Dec2019
%
%             2.  The CSV file CACL_SubjMorphometry_11Dec2019.csv must
%             be in the current directory.
%
%             3.  M-files comp_msh.m, get_subj.m and nod2ele.m must be
%             in the current path or directory.
%
%             4.  For sagittal and coronal slice data.
%
%             5.  To get the MS-Excel spreadsheets, the program must be
%             run from a MS-Windows OS with MS-Excel.
%
%     11-Dec-2019 * Mack Gardner-Morse
%

%#######################################################################
%
% Get Results Subdirectories
%
ddir = fullfile('..','CACL MRI Master');
ddir0 = fullfile(ddir,'TibiaCartThk_21_13Dec2019');
ddir1 = fullfile(ddir,'TibiaCartThkR_21_13Dec2019');
%
% Subject Demographics (Morphometry) CSV File Name
%
dnam = 'CACL_SubjMorphometry_11Dec2019.csv';
%
% Headers for Output File
%
hdrs = {'SID','ID#','Leg','Sex','Age(yrs)','Ht(cm)','Wt(kg)','BMI', ...
        'Comp','Region','Trial'};
nhdrs = size(hdrs,2);   % Number of headers
dcol = char(nhdrs+65);  % Column index to first column of thickness data
%
% Load Grid Data
%
gnam = 'tgrid08_21.mat';
griddat = load(fullfile(ddir0,gnam),'nl','ns','quadl','xql','yql');
nl = griddat.nl;        % Number of points in lateral compartment
ns = griddat.ns;        % Number of subjects
quadl = griddat.quadl;  % Grid quadrilateral connectivity
xql = griddat.xql;      % X coordinates of lateral compartment
yql = griddat.yql;      % Y coordinates of lateral compartment
%
% Load Cartilage Thicknesses
%
cnam = 'tcart08_21.mat';
dat0 = load(fullfile(ddir0,cnam),'cthklc','cthkls','ilegs','kids');
cthklc0 = dat0.cthklc;
cthkls0 = dat0.cthkls;
ilegs = uint8(dat0.ilegs);
%
dat1 = load(fullfile(ddir1,cnam),'cthklc','cthkls','kids');
cthklc1 = dat1.cthklc;
cthkls1 = dat1.cthkls;
%
% Check Subject IDs
%
kids0 = sortrows(dat0.kids);
if ~strcmp(kids0,sortrows(dat1.kids))
  error([' *** ERROR in cthkwri8r:  Not the same subjects at the ', ...
         'two time points!']);
end
%
% Get Subject Numbers
%
[~,ids] = get_subj(kids0);
%
% Get Demographics (Morphometry)
%   Column 1 - Subject identification numbers (IDs)
%   Column 2 - Subject sex numbers (0 - female, 1 - male)
%   Column 3 - Subject ages in years
%   Column 4 - Subject heights in cm
%   Column 5 - Subject weights in kg
%   Column 6 - Subject body mass indexes (BMI) in kg/m^2
%
demog = csvread(dnam,1);
idsall = demog(:,1);
sahwb = demog(:,2:6);
[~,idx] = ismember(ids,idsall);
sahwb = sahwb(idx,:);
%
% Set NaNs to Negative One (-1) in Cartilage Thicknesses
%
idlc0 = isnan(cthklc0);
idls0 = isnan(cthkls0);
idlc1 = isnan(cthklc1);
idls1 = isnan(cthkls1);
%
cthklc0(idlc0) = -1;
cthkls0(idls0) = -1;
cthklc1(idlc1) = -1;
cthkls1(idls1) = -1;
%
% Average Thicknesses
%
idlc0 = ~idlc0;
idls0 = ~idls0;
idlc1 = ~idlc1;
idls1 = ~idls1;
%
cthkl0 = (idlc0.*cthklc0+idls0.*cthkls0)./(idlc0+idls0);
cthkl1 = (idlc1.*cthklc1+idls1.*cthkls1)./(idlc1+idls1);
%
cthkl0(isnan(cthkl0)) = -1;
cthkl1(isnan(cthkl1)) = -1;
%
% Find Thicknesses Within Grid and Compress Grid
%
idc = any(idlc0')'|any(idls0')'|any(idlc1')'|any(idls1')';
%
idcxl = [min(xql(idc)) max(xql(idc))];
idx = xql>=idcxl(1)&xql<=idcxl(2);
%
idcyl = [min(yql(idc)) max(yql(idc))];
idy = yql>=idcyl(1)&yql<=idcyl(2);
%
idc = find(idx&idy);
%
idq = nod2ele(idc,quadl,3);            % Index to compressed rectangle connectivity
%
[xyl,quadlc] = comp_msh([xql yql],quadl(idq,:));
%
xqlc = xyl(:,1);
nlc = size(xqlc,1);
yqlc = xyl(:,2);
%
cthkl0 = cthkl0(idc,:);
cthkl1 = cthkl1(idc,:);
%
% Labels for Points in Lateral Mesh
%
lbls = [repmat('Pt_',nlc,1) strjust(int2str((1:nlc)'),'left')];
lbls = cellstr(lbls)';
%
% Define AP Regions
% Based on the AP Dimensions of the Compressed Lateral Tibia Grid with
% Mean Thicknesses (See Plots at the End of the Program)
%
xpos = -15.5;
xant = 7.5;
idxp = find(xqlc<=xpos);% Posterior region
idxa = find(xqlc>=xant);% Anterior region
idregl = zeros(nlc,1);
idregl(idxp) = -1;
idregl(idxa) = 1;
%
% Save Compressed Grid and AP Regions
%
idot = strfind(gnam,'.');
snam = [gnam(1:idot-1) 'c' gnam(idot:end)];
save(fullfile(ddir1,snam),'idc','idq','idregl','nlc','quadlc', ...
     'xant','xpos','xqlc','yqlc');
%
% Write Out to Three (3) MS-Excel spreadsheets:  tcthk_ant_lat.xlsx,
% tcthk_ctr_lat.xlsx and tcthk_pos_lat.xlsx.
%
kids = cellstr(kids0);
ext_nams = ['pos'; 'ctr'; 'ant'];
cmprt_nams = ['lat'; 'med'];
irow = int2str(ns+2);
%
for k = 1:3             % Regions (1 => Posterior, 2 => Central, 3 => Anterior
   xlsnamr = fullfile(ddir1,['tcthk_' ext_nams(k,:)]);
   for l = 1:1          % Compartments (1 => Lateral, 2 => Medial)
      xlsnam = [xlsnamr '_' cmprt_nams(l,:) '.xlsx'];
      if l==1
        idp = find(idregl==k-2);
      else
        idp = find(idregm==k-2);
      end
%
% Write Out Headers
%
      xlswrite(xlsnam,hdrs,1,'A1');
      xlswrite(xlsnam,lbls(idp),1,[dcol '1']);
%
% Write Out First Trial (0) Data
%
      xlswrite(xlsnam,kids,1,'A2');    % Subject ID
      xlswrite(xlsnam,ids,1,'B2');     % Subject ID number
      xlswrite(xlsnam,ilegs,1,'C2');   % Leg (Left=0/Right=1)
      xlswrite(xlsnam,sahwb,1,'D2');   % Demographics/Morphometry
      xlswrite(xlsnam,repmat(l-1,ns,1),1,'I2');  % Compartment (Lat=0)
      xlswrite(xlsnam,repmat(k-2,ns,1),1,'J2');  % Region (Post=-1/Ctr=0/Ant=1)
      xlswrite(xlsnam,zeros(ns,1),1,'K2');       % Trial
      if l==1
        xlswrite(xlsnam,cthkl0(idp,:)',1,[dcol '2']);% Thicknesses
      else
        xlswrite(xlsnam,cthkm0(idp,:)',1,[dcol '2']);% Thicknesses
      end
%
% Write Out Second Trial (1) Data
%
      xlswrite(xlsnam,kids,1,['A' irow]);   % Subject ID
      xlswrite(xlsnam,ids,1,['B' irow]);    % Subject ID number
      xlswrite(xlsnam,ilegs,1,['C' irow]);  % Left==0/Right==1
      xlswrite(xlsnam,sahwb,1,['D' irow]);  % Demographics/Morphometry
      xlswrite(xlsnam,repmat(l-1,ns,1),1,['I' irow]); % Compartment (Lat=0)
      xlswrite(xlsnam,repmat(k-2,ns,1),1,['J' irow]); % Region (Post=-1/Ctr=0/Ant=1)
      xlswrite(xlsnam,ones(ns,1),1,['K' irow]);       % Trial
      if l==1
        xlswrite(xlsnam,cthkl1(idp,:)',1,[dcol irow]);% Thicknesses
      else
        xlswrite(xlsnam,cthkm1(idp,:)',1,[dcol irow]);% Thicknesses
      end
   end
end
%
% Generate Plots of the Differences of the Means
%
id0 = ~(cthkl0==-1);    % Valid thicknesses
id1 = ~(cthkl1==-1);    % Valid thicknesses
cthkl0m = sum(id0.*cthkl0,2)./sum(id0,2);   % Mean thicknesses
cthkl1m = sum(id1.*cthkl1,2)./sum(id1,2);   % Mean thicknesses
%
cthklmd = cthkl1m-cthkl0m;
%
figure;
orient tall;
patch(xqlc(quadlc'),yqlc(quadlc'),cthklmd(quadlc'), ...
      cthklmd(quadlc'),'FaceColor','interp','EdgeColor','interp');
%
hold on;
id = ~isnan(cthklmd);
plot3(xqlc(id),yqlc(id),10*ones(sum(id),1),'k.','MarkerSize',7);% Grid
plot3([7.5; 7.5],[0 40],[11 11],'r-','LineWidth',2);  % Anterior
plot3([-15.5; -15.5],[0 40],[11 11],'r-','LineWidth',2);   % Posterior
%
view(-90,90);
axis tight;
axis equal;
axis off;
%
colorbar;
caxis([-1.0 1.0]);
colormap(jet);
%
title({'Differences of the Means';'\fontsize{18} Anterior'}, ...
      'FontSize',24,'FontWeight','bold');
%
psnam = 'tcthkwri8r.ps';
print('-dpsc2','-fillpage','-r600',fullfile(ddir1,psnam));
%
% Generate Plots of the Means of the Differences
%
cthkld = cthkl1.*id1-cthkl0.*id0;
idd = id0&id1;
cthkldm = sum(cthkld,2)./sum(idd,2);
idi = isinf(cthkldm);
cthkldm(idi) = NaN;
%
figure;
orient tall;
patch(xqlc(quadlc'),yqlc(quadlc'),cthkldm(quadlc'), ...
      cthkldm(quadlc'),'FaceColor','interp','EdgeColor','interp');
%
hold on;
id = ~isnan(cthkldm);
plot3(xqlc(id),yqlc(id),10*ones(sum(id),1),'k.','MarkerSize',7);% Grid
plot3([7.5; 7.5],[0 40],[11 11],'r-','LineWidth',2);  % Anterior
plot3([-15.5; -15.5],[0 40],[11 11],'r-','LineWidth',2);   % Posterior
%
view(-90,90);
axis tight;
axis equal;
axis off;
%
colorbar;
caxis([-1.0 1.0]);
colormap(jet);
%
title({'Means of the Differences';'\fontsize{18} Anterior'}, ...
      'FontSize',24,'FontWeight','bold');
%
print('-dpsc2','-fillpage','-r600','-append',fullfile(ddir1,psnam));
%
return