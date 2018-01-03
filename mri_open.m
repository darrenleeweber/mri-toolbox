function [mri] = mri_open(mri)

% mri_open - function to call various mri data tools
%
% Useage: [mri] = mri_open(mri)
%
% mri is a parameter structure (see mri_toolbox_defaults for
% more details). In this function, it should contain at least
% the following string fields:
%       
%       mri.path - the directory location of the file to load
%       mri.file - the name of the file to load
%       mri.type - the file format (Analyze, FreeSurfer)
%       
%       Analyze is documented in avw_read etc.
%       FreeSurfer: http://surfer.nmr.mgh.harvard.edu/
%       
% The return structure creates or updates mri.data, which contains:
%       
%       mri.data.hdr     struct, eg see avw_hdr_read
%       mri.data.img     3D matrix of image values
%       
% To plot the data returned, set mri.plot = 1 before loading, or use:
%       
%       avw_view(mri.data)
%       
% See also, avw_img_read, cor_img_read, avw_view
% 

% $Revision: 1.2 $ $Date: 2004/02/07 01:41:51 $

% Licence:  GNU GPL, no express or implied warranties
% History:  08/2002, Darren.Weber@flinders.edu.au
%           11/2002, Darren.Weber@flinders.edu.au
%                    corrected some bugs and mistakes on mri.type
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

version = '$Revision: 1.2 $';
fprintf('MRI_OPEN [v %s]\n',version(11:15));

if ~exist('mri','var'),
    mri = mri_toolbox_defaults;
    fprintf('...creating default mri structure.\n');
elseif isempty(mri),
    mri = mri_toolbox_defaults;
    fprintf('...creating default p structure.\n');
end

[path,name,ext] = fileparts(strcat(mri.path,filesep,mri.file));
file = fullfile(path,[name ext]);

type = lower(mri.type);

switch type,
    
    case 'analyze',
        
        fprintf('...loading Analyze MRI from:\n... %s\n\n',file);
        
        % see avw_img_read for details about orientation
        switch mri.orient
            case 'auto',                mriOrient = '';
            case 'axial unflipped',     mriOrient = 0;
            case 'coronal unflipped',   mriOrient = 1;
            case 'sagittal unflipped',  mriOrient = 2;
            case 'axial flipped',       mriOrient = 3;
            case 'coronal flipped',     mriOrient = 4;
            case 'sagittal flipped',    mriOrient = 5;
            otherwise,                  mriOrient = '';
        end
        
        [ mri.data, mri.IEEEMachine ] = avw_read(file, mriOrient, mri.IEEEMachine);
        
    case 'brainstorm',
        
        fprintf('...BrainStorm not supported yet\n\n');
        return
        %fprintf('...loading BrainStorm data from:\n... %s\n',file);
        
    case {'cor','freesurfer','freesurfer cor'},
        
        % Get Freesurfer data
        [ mri.data, mri.IEEEMachine ] = cor_img_read(path, mri.IEEEMachine);
        
    case 'ge',
        
        % extract series number from path
        separators = findstr(mri.path,filesep);
        seriesPath = mri.path(1:separators(end-1)-1);
        seriesN = mri.path(separators(end-1)+1:separators(end)-1);
                
        % Get GE data
        [ mri.data, mri.IEEEMachine ] = ge_series_read(seriesPath, seriesN);
        if mri.plot,
            fprintf('...cannot plot GE data as yet, use ge_series2avw and avw_view\n');
            mri.plot = 0;
        end
        
    otherwise,
        fprintf('...MRI format: %s\n', mri.type);
        fprintf('...Sorry, cannot load this data format at present.\n\n');
        return;
end

if mri.plot,
    avw_view(mri.data);
end

return
