function [ NeighborList ] = GetNeighborList( Pos, h, Domain_x, Domain_y, NumFluidParticle )
%Function based off of NPP by Siva Srinivas Kolukula, PhD
%Taken from http://www.mathworks.com/matlabcentral/fileexchange/44334-nearest-neighboring-particle-search-using-all-particles-search
%
%   Modified by Tim Voss and Felipe Gutierrez
%
%
%Input:
%   Pos, Position of particles in form [x1, y1; x2, y2;...]
%   h, iteraction radius
%   Domain_x, maximum x coordinate to be meshed
%   Domain_y, maximum y coordinate to be meshed
%
%Output:
%   Neighborlist: an array containing an array of all particles in the
%   vacinity, not organized, includes particle in question

%Create mesh from 0 coordinate to maximum domain by domain/h

%Find Domain X and Y size, remove user input

x = linspace(0,Domain_x,Domain_x/h);
y = linspace(0,Domain_y,Domain_y/h);
[X Y] = meshgrid(x,y);

coor = Pos' ;
List = cell(NumFluidParticle,1) ;                 %Allocate list memory

h = h*2.000001;
for i = 1:NumFluidParticle
    % Get the distance bw ith point and rest all points
    data = repmat(coor(i,:),[length(coor),1])-coor ;
    dist = sqrt(data(:,1).^2+data(:,2).^2);
    % Arrange the distances in ascending order
    [val pos] = sort(dist) ;
    % Pick the points which lie within 'h'
    neighbour = pos(val<=h) ;
    %Output
    List{i} = neighbour ; 
end

NeighborList = List;
end

