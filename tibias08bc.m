%#######################################################################
%
%               * TIBIAS 08 Bone and Cartilage Program *
%
%          Program to run tibias08b.m and tibias08c.m in all selected
%     subject folders.  See tibias08b.m and tibias08c.m for more
%     information.
%
%     NOTES:  1. The Matlab M-files fix_pts.m, li_clos.m, mk_tri6.m,
%             plane_fit.m, rd_roi4.m, rotxyz.m, tibia_cs8.m,
%             tri_area.m, tri_fix2.m and tri_norm.m must be in the
%             current directory or path.
%
%             2.  The input CSV files must follow a specific naming
%             convention.  See tibias08b.m and tibias08c.m.
%
%     25-Nov-2019 * Mack Gardner-Morse
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
%
% Loop through Subjects (Tibias)
%
for ilp = 1:ns
%
% Get Subject Directory
%
   pnam = fullfile(ddir,snams(ilp,:)); % Subject directory
%
% Get Bone File Name
%
   bfnam = fullfile(pnam,'*_SAG_TIB*.csv');      % Bone CSV file names
   d = dir(bfnam);
   if isempty(d)
     error([' *** ERROR in tibias08bc:  Not able to find tibia ', ...
            'bone file:\n\n   ',bfnam,'\n']);
   else
     fnam = d.name;     % Bone CSV file name
   end
%
% Run TIBIAS08B
%
   tibias08b;
%
   clear bfnam fnam;
   close all;
   fclose all;
   clc;
%
% Get Cartilage File Name
%
   cfnam = fullfile(pnam,'*_SAGAR_TIB*.csv');    % Cartilage CSV file names
   d = dir(cfnam);
   if isempty(d)
     error([' *** ERROR in tibias08bc:  Not able to find tibia ', ...
            'cartilage file:\n\n   ',cfnam,'\n']);
   else
     fnams = d.name;    % Cartilage CSV file name
   end
%
% Run TIBIAS08C
%
   tibias08c;
%
   clear cfnam fnams;
   close all;
   fclose all;
   clc;
%
end
%
return