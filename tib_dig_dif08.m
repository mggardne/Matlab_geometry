%#######################################################################
%
%             * TIBia DIGitization DIFFerence Program 08 *
%
%          Program to calculate and plot the Z difference between the
%     coronal and sagittal MRI digitizations of the tibia plateau
%     cartilage.
%
%     NOTES:  1.  The names of the regions of interest (ROIs) must
%             follow a specific convention.
%
%             2.  The M-files dig_zdif3.m, meshbnd4.m, mk_tri4.m,
%             plane_fit.m, rd_roi4.m, sl_info.m and tri_fix2.m must
%             be in the current path or directory.
%
%     02-Jul-2019 * Mack Gardner-Morse
%

%#######################################################################
%
% Compartment Labels
%
tcmpt = ['Lateral'
         'Medial '];
%
% Get Sagittal Cartilage CSV File Name
%
[fnams,pnam] = uigetfile({'*_SAGAR_TIB*.csv', ...
                         'Sagittal tibial cartilage CSV files'}, ...
                         ['Please Select Tibial Sagittal ', ...
                         'Cartilage CSV Files']);
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
  error(' *** ERROR in tib_dig_dif08:  Unable to find coronal file!');
end
%
if size(fnamc,1)>1
  warning([' *** WARNING in tib_dig_dif08:  More than one ', ...
           'coronal file!']);
end
%
fnamc = fnamc(1).name;
%
% Get Tibia Cartilage Data
%
rois = rd_roi4(fullfile(pnam,fnams));  % Sagittal data
roic = rd_roi4(fullfile(pnam,fnamc));  % Coronal data
%
roinams = upper(char(rois.name));
if size(roinams,1)<2
  error(' *** ERROR in tib_dig_dif08:  Missing compartment name!');
end
%
ids(2) = strmatch('MACS',roinams,'exact');  % Medial compartment
ids(1) = strmatch('LACS',roinams,'exact');  % Lateral compartment
%
roinamc = upper(char(roic.name));
if size(roinamc,1)<2
  error(' *** ERROR in tib_dig_dif08:  Missing compartment name!');
end
%
idc(2) = strmatch('MACC',roinamc,'exact');  % Medial compartment
idc(1) = strmatch('LACC',roinamc,'exact');  % Lateral compartment
%
% Setup Differences Figures
%
hf1 = figure;           % Lateral difference figure
hf2 = figure;           % Medial difference figure
%
% Loop through Compartments
%
for l = 1:2
%
% Get Data
%
   dat1 = rois(ids(l)).data';
   dat2 = roic(idc(l)).data';
%
% Get Differences
%
   if l==1
     [davg,dstd,dmin,dmax,daavg,dif] = dig_zdif3(dat1,dat2,true,hf1);
   else
     [davg,dstd,dmin,dmax,daavg,dif] = dig_zdif3(dat1,dat2,true,hf2);
   end
%
end
%
% Update Titles on Differences Figures
%
figure(hf1);
title({[kid ' - Lateral']; 'Digitization Differences'}, ...
      'FontSize',16,'FontWeight','bold','Interpreter','none');
%
figure(hf2);
title({[kid ' - Medial']; 'Digitization Differences'}, ...
      'FontSize',16,'FontWeight','bold','Interpreter','none');
%
return