function [vertices] = mri_scalp

% mri_scalp - Find the scalp surface of an MRI volume
%
% A very course isosurface routine, unlikely to produce
% reliable results.
% 
% See also mesh_shrinkwrap
% 

% $Revision: 1.2 $ $Date: 2004/02/07 01:41:51 $

% Licence:  GNU GPL, no express or implied warranties
% History:  03/2002, Darren.Weber@flinders.edu.au
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load mri;

D = squeeze(D);

Ds = smooth3(D);
FV = isosurface(Ds,5);

vertices = FV.vertices;


% delaunayn can be used to tesselate these vertices
% dsearchn can be used to identify nearest neighbours

r = 10000;
NFV = reducepatch(FV,r); % reduces the faces of struct fv.


Nvertices = size(NFV.vertices,1);

fprintf('\n...Found %d vertices.\n\n',Nvertices);

plot = 1;
if isequal(plot,1),
    P = patch(NFV,'FaceColor',[1,.75,.65],'EdgeColor','none');
    isonormals(Ds,P);
    view(45,30), axis tight, daspect([1,1,.4]), rotate3D;
    lightangle(45,30); 
    set(gcf,'Renderer','zbuffer'); lighting phong
    set(P,'SpecularColorReflectance',0,'SpecularExponent',50);
end


return
