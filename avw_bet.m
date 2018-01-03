function [ u ] = avw_bet(avw,bt)

% avw_bet
%
% This code is developed entirely on the basis of Smith, S. (2002). Fast robust automated brain extraction. Human Brain Mapping, 17(3): 143-155.



% see Smith (2002, p. 150)
% bt can vary between 0-1
if ~exist('bt','var'), bt = 0.5; end
if isempty(bt), bt = 0.5; end

fprintf('...estimating brain surface...'); tic;

[FV,t,t02,t98,tm,COG] = avw_bet_init(avw);

% this loops 1000 times

% component 1
[S,Sn,St] = avw_bet_update1(FV);

% component 2
[f2,sn,sn_unit,L] = avw_bet_update2(FV,COG.voxels);

% component 3
%[ f3, t1, Imin, Imax] = avw_bet_update3(avw,bt,v,sn)
f3 = avw_bet_update3(avw,bt,v,sn,t,t02,tm);

% Nvert is size(FV.vertices,1)
% st is a vector Nvert x 3 (tangential to local surface)
% sn is a vector Nvert x 3 (normal to local surface)
% sn_unit is the unit vector of sn, a vector Nvert x 3
% f3 is a scalar Nvert x 1 (fractional update of intensity)
% L  is a scalar Nvert x 1 (the mean distance of a vertex to its neighbours)

% check on the units of sn and L - voxels or mm?

u = (0.5 * st) + (f2 .* sn) + (0.05 .* f3 .* L .* sn_unit);


time = toc; fprintf('done (%5.2f sec)\n',time);
return






%------------------------------------------------------
function [FV,t,t02,t98,tm,COG] = avw_bet_init(avw)

% t02 and t98 are 2% and 98% of cumulative volume histogram
[bins,freq,freq_nozero] = avw_histogram(avw,2);

% this seems to work, but not sure if it is a correct definition of
% 'cumulative volume histogram'
binbyfreq = bins .* freq;

t02_tmp = 0.02 * sum(binbyfreq);
t98_tmp = 0.98 * sum(binbyfreq);

i = 1;
while sum(binbyfreq(1:i)) < t02_tmp, i = i + 1; end
t02 = bins(i-1);

i = 1; while sum(binbyfreq(1:i)) < t98_tmp, i = i + 1; end
t98 = bins(i);

% t is the 'brain/background' threshold, which is used to roughly estimate
% the position of the center of gravity of the brain/head in the image
% volume

% t is simply set to lie 10% of the way between t02 and t98
t = (0.1 * (t98 - t02)) + t02;

%Create binarized image to determine Center of gravity (COG) and radius of
%sphere volume to determine tm.
bin = avw;
% find all voxels with intensity less than t
voxels_lt_t = find(bin.img < t);
bin.img(voxels_lt_t) = 0;
% find all voxels with intensity greater than t98
voxels_gt_t98 = find(bin.img > t98);
bin.img(voxels_gt_t98) = 0;
% binary image
nonzero = find(bin.img);
bin.img(nonzero) = 1;

% find weighted sum of positions (CENTER OF GRAVITY (COG))
COG = avw_center_mass(bin);
xCOG = COG.voxels(1);
yCOG = COG.voxels(2);
zCOG = COG.voxels(3);

% now find the radius for calculation of tm and calculate the radius of all
% voxel locations with respect to COG to identify all the voxels within the
% radius of the sphere used to calculate tm.

%this might be useful with sperical source modeling.

%1. loop over voxels in bin.img...over thresholded volume or over image?
%2. calculate radius from Center for each radius
%3. the mean of the largest 10% of these radii = estimated radius


% double check this
xi = 1:avw.hdr.dime.dim(2);
yj = 1:avw.hdr.dime.dim(3);
zk = 1:avw.hdr.dime.dim(4);
voxel_radius(i,j,k)= sqrt( (xi-xCOG).^2 + (yj-yCOG).^2 + (zk-zCOG).^2 );

max_voxel_radius = max(max(max(voxel_radius)));
thresh_voxel_radius = 0.9 * max_voxel_radius;
voxel_radius_gt_90_index = find(voxel_radius > thresh_voxel_radius);
mean_radius = mean(voxel_radius_gt_90(voxel_radius_gt_90_index));

% "Finally, the median intensity of all points within a sphere of the
% estimated radius and centered on the estimated COG is found (tm)"
% Smith(2002, p. 146).

index_voxels_inside_mean_radius = find(voxel_radius <= mean_radius);
tm = median(avw.img(index_voxels_inside_mean_radius));

% Initialise the spherical tesselation
% with 4 recursions we get 2562 vertices,
% with 5 recursions we get 10,242 vertices;
FV = sphere_tri('ico',4,mean_radius);
FV.vertices = FV.vertices + [xCOG,yCOG,zCOG];

return







%------------------------------------------------------
function [S,Sn,St] = avw_bet_update1(FV),

% This function implements Smith, S. (2002), Fast robust automated brain
% extraction.  Human Brain Mapping, 17: 143-155.  It corresponds to update
% component 1: within surface vertex spacing.

if isfield(FV,'edge'),
    if isempty(FV.edge),
        FV.edge = mesh_edges(FV);
    end
else
    FV.edge = mesh_edges(FV);
end

[normals,unit_normals] = mesh_vertex_normals(FV);

Nvert = size(FV.vertices,1);

for index = 1:Nvert,

    v = FV.vertices(index,:);
    x = FV.vertices(index,1);
    y = FV.vertices(index,2);
    z = FV.vertices(index,3);

    unit_normal = unit_normals(index,:);

    % Find neighbouring vertex coordinates
    vi = find(FV.edge(index,:));  % the indices
    neighbour_vertices = FV.vertices(vi,:);
    X = neighbour_vertices(:,1);
    Y = neighbour_vertices(:,2);
    Z = neighbour_vertices(:,3);

    % Find neighbour mean location; this is 'mean position of A and B' in
    % figure 4 of Smith (2002)
    Xmean = mean(X);
    Ymean = mean(Y);
    Zmean = mean(Z);

    % Find difference in distance between the vertex of interest and its
    % neighbours; this value is 's' and 'sn' in figure 4 of
    % Smith (2002, eq. 1 to 4)
    s = [ Xmean - x, Ymean - y, Zmean - z]; % inward toward mean

    % Find the vector sn
    sn = dot( s, unit_normal ) * unit_normal;

    % Find the vector st
    st =  s - sn; % absolute value

    S(index,:) = s;
    Sn(index,:) = sn;
    St(index,:) = st;

end

return






%------------------------------------------------------
function [f2,sn,sn_unit,L] = avw_bet_update2(FV,origin),

% This function adapts Smith, S. (2002), Fast robust automated brain
% extraction.  Human Brain Mapping, 17: 143-155.  This function
% corresponds to update component 2: surface smoothness control.

xo = origin(1); yo = origin(2); zo = origin(3);

Nvert = size(FV,1);

f2 = zeros(Nvert,1);
sn = zeros(Nvert,3);

for index = 1:Nvert,

    v = FV.vertices(index,:);
    x = FV.vertices(index,1);
    y = FV.vertices(index,2);
    z = FV.vertices(index,3);

    % Find radial distance of vertex from origin
    r = sqrt( (x-xo)^2 + (y-yo)^2 + (z-zo)^2 );

    % Calculate unit vector
    v_unit_vector = ( v - origin ) / r;

    % Find direction cosines for line from center to vertex
    l = (x-xo)/r; % cos alpha
    m = (y-yo)/r; % cos beta
    n = (z-zo)/r; % cos gamma

    % Find neighbouring vertex coordinates
    vi = find(FV.edge(index,:));  % the indices
    neighbour_vertices = FV.vertices(vi,:);
    X = neighbour_vertices(:,1);
    Y = neighbour_vertices(:,2);
    Z = neighbour_vertices(:,3);

    % Find neighbour radial distances
    r_neighbours = sqrt( (X-xo).^2 + (Y-yo).^2 + (Z-zo).^2 );
    r_neighbours_mean = mean(r_neighbours);
    L(index) = r_neighbours_mean;
    
    % Find difference in radial distance between the vertex of interest and its
    % neighbours; this value approximates the magnitude of sn in
    % Smith (2002, eq. 1 to 4)
    r_diff = r - r_neighbours_mean;

    % Find the vector sn, in the direction of the vertex of interest, given the
    % difference in radial distance between vertex and mean of neighbours
    sn(index,:) = r_diff * v_unit_vector;
    
    snx = sn(index,1);
    sny = sn(index,2);
    snz = sn(index,3);
    
    sn_unit(index,:) = sn(index,:) ./ sqrt( (snx - v(1)).^2 + (sny - v(2)).^2 + (snz - v(3)).^2 );

    % Find distances between vertex and neighbours, using edge lengths.
    % The mean value is l in Smith (2002, eq. 4)
    edge_distance = FV.edge(index,vi);
    edge_distance_mean = mean(edge_distance);

    % Calculate radius of local curvature, solve Smith (2002, eq. 4)
    if r_diff,
        radius_of_curvature = (edge_distance_mean ^ 2) / (2 * r_diff);
    else
        radius_of_curvature = 10000;
    end

    % Define limits for radius of curvature
    radius_min =  3.33; % mm
    radius_max = 10.00; % mm

    % Sigmoid function parameters,
    % "where E and F control the scale and offset of the sigmoid"
    E = mean([(1 / radius_min),  (1 / radius_max)]);
    F = 6 * ( (1 / radius_min) - (1 / radius_max) );

    f2(index) = (1 + tanh( F * (1 / radius_of_curvature - E))) / 2;

%     % multiply sigmoid function by sn
%     move_vector(index,:) = f2 * sn(index,:);
% 
%     FV.vertices(index,:) = v + move_vector;
end

return






% --------------------------------------------------------
function [ f3, t1, Imin, Imax] = avw_bet_update3(avw,bt,v,sn,t,t02,tm)
% component 3 of Smith (2002)


%along a line pointing inward from the current vertex, min and max
%intensities are found. d1 determines how far into the brain the minimum
%intensity is searched for. d2 determines the same for the max intensity.

d1 = 20; % 20 mm (Smith)
d2 = d1/2;

% get intensities along a line of length d1 parallel with the vertex normal
[normals,unit_normals] = mesh_vertex_normals(FV);

Nvert = size(FV.vertices,1);

t1 = zeros(Nvert,1);
f3 = zeros(Nvert,1);
Imin = zeros(Nvert,1);
Imax = zeros(Nvert,1);

for v = Nvert,

    x = FV.vertices(v,1);
    y = FV.vertices(v,2);
    z = FV.vertices(v,3);

    % we need to find points along a line from the vertex point, parallel
    % with the vertex normal

    % vertex location (P) and unit vector of surface normal (n_hat), then
    % points (S) along line parallel to surface normal are obtained with vector
    % addition, such that S = P + (x * n_hat)
    d1_distances = -d1:0.5:0;
    d2_distances =   0:0.5:d2;

    % check the nature of this multiplication (do we need a for loop?)
    S1 = [x y z] + (d1_distances .* unit_normals(v,:));
    S2 = [x y z] + (d2_distances .* unit_normals(v,:));

    X1 = S1(:,1);
    Y1 = S1(:,2);
    Z1 = S1(:,3);

    X2 = S2(:,1);
    Y2 = S2(:,2);
    Z2 = S2(:,3);

    %intensity_d1 = set of intensities at a discrete number of points along
    %this line : [I(0),I(1),I(2),...,I(d1)]
    %intensity_d2 = similar to intensity_d1

    % interpolate volume values at these points
    % ( not sure why have to swap XI,YI here )
    intensity_d1(v,:) = interp3(avw.img,Y1,X1,Z1,'*nearest');
    intensity_d2(v,:) = interp3(avw.img,Y2,X2,Z2,'*nearest');

    Imin(v) = max( t02, min(tm,intensity_d1) );
    Imax(v) = min(  tm, max(t, intensity_d2) );

    % for certain image intenisty distributions it can be varied (in the range
    % 0-1 to give optimal results, but the necessity for this is rare.

    t1(v) = (Imax(v) - t02) * bt + t02;

    f3(v) = ( 2 * (Imin(v) - t1(v)) ) / ( Imax(v) - t02 );

end

return
