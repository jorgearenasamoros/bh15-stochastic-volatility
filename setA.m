function A = setA(vecA)
% SETA  Contemporaneous matrix of the bivariate labor-market model.
% From the Baumeister and Hamilton (2015) replication files.
%   vecA = [beta; alpha],  A = [-beta 1; -alpha 1]
A = [-vecA(1) 1; -vecA(2) 1];
end
