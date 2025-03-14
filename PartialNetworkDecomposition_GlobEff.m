function [PND, M] = PartialNetworkDecomposition_GlobEff(netx, nety)

%%PARTIALNETWORKDECOMPOSITION Decomposes global efficiency between two networks
%
%    PND = PARTIALNETWORKDECOMPOSITION(NETX, NETY) computes the partial network
%    decomposition of global efficiency (i.e. inverse shortest path length) for
%    two binary undirected networks NETX and NETY. Returns a structure with the
%    redundant, unique, and synergistic contributions corresponding to both
%    networks.
%
%    [PND, M] = PARTIALNETWORKDECOMPOSITION(NETX, NETY) returns a matrix M
%    representing the dominant 'mode' of each pair of nodes, where M(i,j) = -1
%    for disconnected pairs, 0 for redundant pairs, 1 (resp.  2) for pairs with
%    unique information from network X (resp. Y), and 3 for synergistic pairs.
%
% This function includes code from the Brain Connectivity Toolbox by Rubinov et
% al: https://sites.google.com/site/bctnet/
%
% Pedro Mediano and Andrea Luppi, January 2023

%% Parameter checks
if ~(ismatrix(netx) && ismatrix(nety) && ...
    size(netx,1) == size(netx,2) && all(size(netx) == size(nety)))
  error("Networks must be square matrices of the same size");
end
if ~(issymmetric(netx) && issymmetric(nety))
  error("Both matrices must be symmetric (i.e. undirected networks).");
end
isxbinary = islogical(netx) || length(unique(netx(:))) <= 2;
isybinary = islogical(nety) || length(unique(nety(:))) <= 2;
if ~(isxbinary && isybinary)
  error("Both inputs must be binary networks (i.e. false/true or 0/1 edges).");
end

% Binarise networks and convert to double
netx = double(netx ~= 0);
nety = double(nety ~= 0);


%% Compute efficiencies in both single networks and their union
Ex  = 1./distance_bin(netx);
Ey  = 1./distance_bin(nety);
Exy = 1./distance_bin(netx | nety);

% Turn diagonal entries to NaN's
D = length(netx);
Ex(1:D+1:end)  = nan;
Ey(1:D+1:end)  = nan;
Exy(1:D+1:end) = nan;


%% Compare these efficiencies to compute PND atoms
PND = [];
PND.red = min(Ex, Ey);
PND.unx = Ex - PND.red;
PND.uny = Ey - PND.red;
PND.syn = Exy - max(Ex, Ey);

% Check that the atoms sum up to the right things
assert(isequalntol( Ex, PND.red + PND.unx));
assert(isequalntol( Ey, PND.red + PND.uny));
assert(isequalntol(Exy, PND.red + PND.unx + PND.uny + PND.syn));


%% If requested, compute the dominant mode of each path
if nargout > 1
  M = -1*ones(size(netx));
  M(PND.syn > 0) = 3;
  M(PND.uny > 0 & M < 0) = 2;
  M(PND.unx > 0 & M < 0) = 1;
  M(PND.red > 0 & M < 0) = 0;
  M(1:length(M)+1:end) = nan;
  
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Small utility function to test whether two arrays containing NaNs are equal
% within a certain tolerance level
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = isequalntol(A, B, tol)
  if nargin < 3 || isempty(tol)
    tol = 1e-12;
  end
  nan_equals = isequal(isnan(A), isnan(B));
  if nan_equals
    num_equals = all(abs(A(~isnan(A)) - B(~isnan(B))) < tol, 'all');
  else
    num_equals = false;
  end
  tf = nan_equals && num_equals;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function taken verbatim from the Brain Connectivity Toolbox on Jan 20th 2023
%
%                   https://sites.google.com/site/bctnet/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D=distance_bin(A)
%DISTANCE_BIN       Distance matrix
%
%   D = distance_bin(A);
%
%   The distance matrix contains lengths of shortest paths between all
%   pairs of nodes. An entry (u,v) represents the length of shortest path 
%   from node u to node v. The average shortest path length is the 
%   characteristic path length of the network.
%
%   Input:      A,      binary directed/undirected connection matrix
%
%   Output:     D,      distance matrix
%
%   Notes: 
%       Lengths between disconnected nodes are set to Inf.
%       Lengths on the main diagonal are set to 0.
%
%   Algorithm: Algebraic shortest paths.
%
%
%   Mika Rubinov, U Cambridge
%   Jonathan Clayden, UCL
%   2007-2013

% Modification history:
% 2007: Original (MR)
% 2013: Bug fix, enforce zero distance for self-connections (JC)

A=double(A~=0);                 %binarize and convert to double format

l=1;                            %path length
Lpath=A;                        %matrix of paths l
D=A;                            %distance matrix

Idx=true;
while any(Idx(:))
    l=l+1;
    Lpath=Lpath*A;
    Idx=(Lpath~=0)&(D==0);
    D(Idx)=l;
end

D(~D)=inf;                      %assign inf to disconnected nodes
D(1:length(A)+1:end)=0;         %clear diagonal
end


