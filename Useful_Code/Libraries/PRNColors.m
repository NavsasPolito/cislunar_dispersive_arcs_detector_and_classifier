function [RGBtriplet] = PRNColors(svID)
% This function assigns a color to the different SV. To do this, it gives a
% random RGB triplet to each of the PRN

%--- Set up seed
rng(mod(keyHash(svID),2^32));
%--- Generate random triplet
RGBtriplet = [rand, rand, rand];