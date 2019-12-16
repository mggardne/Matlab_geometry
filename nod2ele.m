function [ir,n] = nod2ele(nodlst,econn,ncn);
%NOD2ELE  Finds the elements connected to a vector of nodes.
%
%         IR = NOD2ELE(NODLST,ECONN) returns the row index to elements
%         connected to the nodes in the vector, NODLST, that are in the
%         element connectivity matrix, ECONN, where each row represents
%         one element.
%
%         [IR,N] = NOD2ELE(NODLST,ECONN) returns the number of elements,
%         N, connected to the nodes in NODLST.
%
%         IR = NOD2ELE(NODLST,ECONN,NCN) returns the row index to the
%         elements connected to the nodes in NODLST with greater
%         than NCN nodes from NODLST in the element.  NCN must be
%         between zero (0) and one less than the number of nodes in
%         each element.  The default is zero (0) (at least one node
%         from the element is in the input list).
%
%         NOTES:  1.  Note that the number of connected nodes from
%                 NODLST must be greater than the number of connected
%                 nodes, NCN.
%
%         07-Nov-2014 * Mack Gardner-Morse
%

%#######################################################################
%
% Check Inputs
%
if (nargin<2)
  error(' *** ERROR in NOD2ELE:  Not enough input arguments!');
end
%
if (nargin<3)
  ncn = 0;
end
%
nne = size(econn,2);    % Number of nodes in each element
if ncn<0|ncn>nne-1
  error([' *** ERROR in NOD2ELE:  Number of connected nodes, NCN,' ...
         ' must be between zero (0) and one less than the number' ...
         ' of nodes in the elements!']);
end
%
% Get Sizes of Input Data
%
nodlst = nodlst(:);
nn = size(nodlst,1);
%
ne = size(econn,1);
%
% Find Connected Elements
%
nnc = false(ne,nne);
for k = 1:nn
   nnc = nnc|econn==nodlst(k);
end
ir = find(sum(nnc,2)>ncn);
n = size(ir,1);
%
return