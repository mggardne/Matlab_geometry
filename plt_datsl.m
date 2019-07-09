function p = plt_datsl(dat,cmstr,lw,colr,lbl)
% PLT_DATSL Plots ROI data one slice at a time in 3D.
%
%          P = PLT_DATSL(DAT) given a column cell array containing
%          three (3) column matrices with slice coordinate point data,
%          DAT, returns a 3D plot with a default black line for each
%          slice.
%
%          P = PLT_DATSL(DAT,CMSTR) given a cell array containing three
%          (3) column matrices with slice coordinate point data in the
%          cell array, DAT, returns a 3D plot with a marker and color
%          designation of CMSTR.
%
%          P = PLT_DATSL(DAT,CMSTR,LW) adds a designated line width to
%          the integer, LW.
%
%          P = PLT_DATSL(DAT,CMSTR,LW,COLR) adds a designated RGB color
%          to the 1x3 vector, COLR.
%
%          25 July 2012 Daniel R. Sturnick
%

%#######################################################################
%
% Check Input Arguments
%
if nargin<2||isempty(cmstr)
    cmstr = 'k';
end
%
if nargin<3||isempty(lw)
    lw = .5;
    cl = false;
end
%
if nargin<4||isempty(colr)
    cl = false;
else
    cl = true;
end
%
if nargin<5||isempty(lbl)
  lbl = false;
end
%
% Check and Plot Data
%
if iscell(dat)
  nslice = length(dat);
  if lbl
    p = zeros(2*nslice,1);
  else
    p = zeros(nslice,1);
  end
%
  for s = 1:nslice
     xyz = dat{s};
     if cl
         l = plot3(xyz(:,1),xyz(:,2),xyz(:,3),cmstr,'LineWidth',lw, ...
             'MarkerSize',10,'Color',colr);
     else
         l = plot3(xyz(:,1),xyz(:,2),xyz(:,3),cmstr,'LineWidth',lw, ...
             'MarkerSize',10);
     end
%
     p = [p; l];
%
     if lbl
       t = text(xyz(1,1),xyz(1,2),xyz(1,3),int2str(s),'FontSize', ...
                12,'FontWeight','bold');
%                 12,'FontWeight','bold','Color','r');
       p = [p; t];
     end
%
     if s==1
       hold on;
     end
%
  end
%
else
  fprintf(1,'\n No slices to plot\n\n');
  p = NaN;
end
%
return