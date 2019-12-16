function [idnum,nums] = get_subj(strs,delim);
%GET_SUBJ Gets the subject number from the start of the rows in a
%         character array.
%
%         IDNUM = GET_SUBJ(STRS) given a character array STRS, gets the
%         subject number from the start of the rows to the first
%         underscore ('_') as characters in the character array IDNUM.
%
%         IDNUM = GET_SUBJ(STRS,DELIM) uses the characters in DELIM to
%         determine the limits of the subject number from the start of
%         the row.
%
%         [IDNUM,NUMS] = GET_SUBJ(STRS) returns a double vector NUMS
%         containing the subject numbers.
%
%         NOTES:  1.  The delimiter DELIM is case sensitive.
%
%        01-Aug-2013 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)||isempty(delim)
  delim = '_';
end
%
if (nargin<1)
  error(' *** ERROR in GET_SUBJ:  Not enough input data!');
end
%
% Loop through Rows Looking for Subject Numbers
%
nrow = size(strs,1);
idnum = cell(nrow,1);
delim = delim(:)';
%
for k = 1:nrow
   id = findstr(strs(k,:),delim);
   if ~isempty(id)
     idnum{k} = strs(k,1:id(1)-1);
   else
     idnum{k} = [];
   end
end
%
% Convert to a Character Array
%
idnum = char(idnum);
%
% Convert to Numbers
%
if nargout>1
  nums = zeros(nrow,1);
  for k = 1:nrow
     nums(k) = eval(idnum(k,:));
  end
end
%
return