function [ind]=fromPairToIndex(p, q)
m=[[1 6 5];
   [6 2 4];
   [5 4 3]];
ind=m(p,q);