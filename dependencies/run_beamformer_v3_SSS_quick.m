function [bf_outs] = run_beamformer_v3_SSS_quick(fwd_method,sourcepos,S,do_bf,max_or_min,MFC, do_SSS)
%% Function to run whichever beamformer you choose
% [bf_outs] = run_beamformer(fwd_method,sourcepos,S,do_bf,max_or_min,MFC)
% [INPUTS]
% fwd_method      = 'single' (single sphere method).
%                 = 'local' (local spheres method).
%                 = 'BEM' (FieldTrip BEM method requires OpenMEEG)
%                 = 'BEM2' (Stenroos BEM - much faster than FiledTrip)
%                 = 'shell' (single shell method).
% sourcepos       = dipole locations (M x 3).
% S               = parameter structure.
% do_bf           = set to 1 if you want to beamform, or 0 if you only want
%                   the leadfields.
% max_or_min      = look at the maximum or minimum T-Stat by setting 'Max'
%                   or 'Min' (or N/A if not used).
% MFC             = set to 1 if mean field correction is used.
%
% do_SSS          = set to 1 if SSS is used
%
% For all models
% S.mri_file      = .mri file location.
% S.meshes_file   = location of meshes
% S.sensor_info   = sensor information structure (.pos & .ors both
%                   Nchansx3).
% S.M             = The mean field correction matrix (used to multiply the
%                   lead fields).
% S.Lin           = Number of internal components used in SSS
% S.Lout          = Number of external components used in SSS
%
% For Local Spheres
% S.skanect_file  = Coregistered file to save local spheres to.
% S.rad_ors       = Radial orientation associated with each position. This
%                   should be Nchansx3. No matter the channel axis, use the
%                   Z orientation.
%
% For BEM
% S.skanect_file  = Coregistered file to save volume and lead fields to.
%
% For beamforming (do all in Tesla)
% S.C             = Data covariance.
% S.Cinv          = Inverse of data covariance.
% S.Noise_Cr      = Noise covariance.
% S.Ca            = Active covariance. (Optional)
% S.Cc            = Control covariance. (Optional)
%
% [OUTPUTS]
% bf_outs.LF      = leadfield matrix (Nchans x 3 x M). In XYZ format (i.e. a
%                   dipole at [1 0 0], [0 1 0], and [0 0 1]).
% bf_outs.Weights = Weights for each source location (Nchans x M).
% bf_outs.T_stat  = T stat for each source location.
% bf_outs.Z_stat  = Z stat for each source location.

[path.mri, filename.mri] = fileparts(S.mri_file);
path.mri = [path.mri filesep];
try
    mri = ft_read_mri([path.mri filename.mri '.mri']);
catch
    try
        mri = ft_read_mri([path.mri filename.mri '.nii']);
    catch
        error('Check MRI file is  in .mri or .nii format or that FieldTrip is added')
    end
end

if ~exist([S.meshes_file]) || ~isfield(S,'meshes_file')
    [mesh_path,~,~] = fileparts(S.meshes_file);
    if ~exist(mesh_path)
        mkdir(mesh_path)
    end
    if ~exist([mesh_path '/segmentedmri.mat'])
        disp('Segmenting MRI...')
        cfg           = [];
        cfg.output    = {'brain','skull','scalp'};
        segmentedmri  = ft_volumesegment(cfg, mri);
        save([mesh_path '/segmentedmri.mat'],'segmentedmri')
        disp('Done!')
    else
        disp('Loading Segmented MRI...')
        load([mesh_path '/segmentedmri.mat'])
        disp('Done!')
    end
    cfg = [];
    cfg.tissue = {'brain','skull','scalp'};
    cfg.numvertices = [3000 2500 2500];
    mesh1 = ft_prepare_mesh(cfg,segmentedmri);
    for n = 1:size(mesh1,2)
        meshes(n).pnt = mesh1(n).pos;
        meshes(n).tri = mesh1(n).tri;
        meshes(n).unit = mesh1(n).unit;
        meshes(n).name = cfg.tissue{n};
        meshes(n) = ft_convert_units(meshes(n),'m');
    end
    save(S.meshes_file,'meshes','segmentedmri')
else
    load(S.meshes_file)
end
origin = mean(meshes(1).pnt);

%%%%%%%%%%% Single Sphere %%%%%%%%%%%
if strcmpi(fwd_method,'single')
    sensor_info.pos = S.sensor_info.pos;
    sensor_info.ors = S.sensor_info.ors;
    sensor_info.label = strsplit(num2str(1:size(sensor_info.pos)));

    [LF_single_shell,L_reshaped] = Forward_OPM_single_FT_inscript(S.meshes_file,sensor_info,sourcepos);
    disp(['size of L_reshaped is ' num2str(size(L_reshaped))]);
    LF1 = L_reshaped.*1e-9;
    inside_idx = find(LF_single_shell.inside);
    LF = zeros(size(LF1));
    for n = size(sourcepos,1):-1:1
        if LF_single_shell.inside(n)
            nn = find(inside_idx == n);
            LF(:,:,n) = LF1(:,:,nn);
        else
            LF(:,:,n) = zeros(size(LF1(:,:,1)));
        end

    end
    clear LF1 L_reshaped
    if MFC % Account for mean field correction
        for n = 1:size(LF,3)
            LF(:,:,n) = S.M*LF(:,:,n);
        end
    end

    if do_SSS % Apply iterative SSS to lead fields
        ft_progress('init', 'text', 'Performing SSS correction...')      % ascii progress bar
        for n = 1:size(LF, 3)
            ft_progress(n/floor(size(LF,3)), 'Processing voxel %d of %d', n, floor(size(LF,3)))
            [sss_outs]=SSS_opm_lf(LF(:,:,n),sensor_info.pos,sensor_info.ors,S.Lin,S.Lout,S.sss_ins);
            LF(:,:,n) = sss_outs.phi_in;
        end
        ft_progress('close')
    end

    for n = size(sourcepos,1):-1:1
        % project to two tangential fields
        R = sourcepos(n,:) - origin;
        [etheta,ephi] = calctangent_all_bf(R);
        tan_ors = [etheta;ephi];
        ltan(:,:,n) = LF(:,:,n)*tan_ors';
    end

    %%%%%%%%%%% Local Spheres %%%%%%%%%%%
elseif strcmpi(fwd_method,'local')
    sensor_info.pos = S.sensor_info.pos;
    sensor_info.ors = S.sensor_info.ors;
    sensor_info.label = strsplit(num2str(1:size(sensor_info.pos)));

    [LF_single_shell,L_reshaped] = Forward_OPM_local_FT_inscript(S.meshes_file,sensor_info,sourcepos);
    disp(['size of L_reshaped is ' num2str(size(L_reshaped))]);
    LF1 = L_reshaped.*1e-9;
    inside_idx = find(LF_single_shell.inside);
    LF = zeros(size(LF1));
    for n = size(sourcepos,1):-1:1
        if LF_single_shell.inside(n)
            nn = find(inside_idx == n);
            LF(:,:,n) = LF1(:,:,nn);
        else
            LF(:,:,n) = zeros(size(LF1(:,:,1)));
        end

    end
    clear LF1 L_reshaped
    if MFC % Account for mean field correction
        for n = 1:size(LF,3)
            LF(:,:,n) = S.M*LF(:,:,n);
        end
    end

    if do_SSS % Apply iterative SSS to lead fields
        ft_progress('init', 'text', 'Performing SSS correction...')      % ascii progress bar
        for n = 1:size(LF, 3)
            ft_progress(n/floor(size(LF,3)), 'Processing voxel %d of %d', n, floor(size(LF,3)))
            [sss_outs]=SSS_opm_lf(LF(:,:,n),sensor_info.pos,sensor_info.ors,S.Lin,S.Lout,S.sss_ins);
            LF(:,:,n) = sss_outs.phi_in;
        end
        ft_progress('close')
    end

    for n = size(sourcepos,1):-1:1
        % project to two tangential fields
        R = sourcepos(n,:) - origin;
        [etheta,ephi] = calctangent_all_bf(R);
        tan_ors = [etheta;ephi];
        ltan(:,:,n) = LF(:,:,n)*tan_ors';
    end

elseif strcmpi(fwd_method,'BEM')
    sensor_info.pos = S.sensor_info.pos;
    sensor_info.ors = S.sensor_info.ors;
    sensor_info.label = strsplit(num2str(1:size(sensor_info.pos)));

    [LF_single_shell,L_reshaped] = Forward_OPM_BEM_FT_inscript(S.meshes_file,sensor_info,sourcepos);
    disp(['size of L_reshaped is ' num2str(size(L_reshaped))]);
    LF1 = L_reshaped.*1e-9;
    inside_idx = find(LF_single_shell.inside);
    LF = zeros(size(LF1));
    for n = size(sourcepos,1):-1:1
        if LF_single_shell.inside(n)
            nn = find(inside_idx == n);
            LF(:,:,n) = LF1(:,:,nn);
        else
            LF(:,:,n) = zeros(size(LF1(:,:,1)));
        end

    end
    clear LF1 L_reshaped
    if MFC % Account for mean field correction
        for n = 1:size(LF,3)
            LF(:,:,n) = S.M*LF(:,:,n);
        end
    end

    if do_SSS % Apply iterative SSS to lead fields
        ft_progress('init', 'text', 'Performing SSS correction...')      % ascii progress bar
        for n = 1:size(LF, 3)
            ft_progress(n/floor(size(LF,3)), 'Processing voxel %d of %d', n, floor(size(LF,3)))
            [sss_outs]=SSS_opm_lf(LF(:,:,n),sensor_info.pos,sensor_info.ors,S.Lin,S.Lout,S.sss_ins);
            LF(:,:,n) = sss_outs.phi_in;
        end
        ft_progress('close')
    end

    for n = size(sourcepos,1):-1:1
        % project to two tangential fields
        R = sourcepos(n,:) - origin;
        [etheta,ephi] = calctangent_all_bf(R);
        tan_ors = [etheta;ephi];
        ltan(:,:,n) = LF(:,:,n)*tan_ors';
    end

elseif strcmpi(fwd_method,'BEM2')
    sensor_info.pos = S.sensor_info.pos;
    sensor_info.ors = S.sensor_info.ors;
    sensor_info.label = strsplit(num2str(1:size(sensor_info.pos)));

    [LF_BEM2,L_reshaped] = Forward_OPM_BEM_stenroos_inscript(S.meshes_file,sensor_info,sourcepos);
    disp(['size of L_reshaped is ' num2str(size(L_reshaped))]);
    LF1 = L_reshaped.*1e-9;
    LF = zeros(size(LF1));
    for n = size(sourcepos,1):-1:1
        LF(:,:,n) = LF1(:,:,n);
    end
    clear LF1 L_reshaped
    if MFC % Account for mean field correction
        for n = 1:size(LF,3)
            LF(:,:,n) = S.M*LF(:,:,n);
        end
    end

    if do_SSS % Apply iterative SSS to lead fields
        ft_progress('init', 'text', 'Performing SSS correction...')      % ascii progress bar
        for n = 1:size(LF, 3)
            ft_progress(n/floor(size(LF,3)), 'Processing voxel %d of %d', n, floor(size(LF,3)))
            [sss_outs]=SSS_opm_lf(LF(:,:,n),sensor_info.pos,sensor_info.ors,S.Lin,S.Lout,S.sss_ins);
            LF(:,:,n) = sss_outs.phi_in;
        end
        ft_progress('close')
    end

    for n = size(sourcepos,1):-1:1
        % project to two tangential fields
        R = sourcepos(n,:) - origin;
        [etheta,ephi] = calctangent_all_bf(R);
        tan_ors = [etheta;ephi];
        ltan(:,:,n) = LF(:,:,n)*tan_ors';
    end

elseif strcmpi(fwd_method,'shell')
    sensor_info.pos = S.sensor_info.pos;
    sensor_info.ors = S.sensor_info.ors;
    sensor_info.label = strsplit(num2str(1:size(sensor_info.pos)));

    [LF_single_shell,L_reshaped] = Forward_OPM_shell_FT_inscript(S.meshes_file,sensor_info,sourcepos);
    disp(['size of L_reshaped is ' num2str(size(L_reshaped))]);
    LF1 = L_reshaped.*1e-9;
    inside_idx = find(LF_single_shell.inside);
    LF = zeros(size(LF1));
    for n = size(sourcepos,1):-1:1
        if LF_single_shell.inside(n)
            nn = find(inside_idx == n);
            LF(:,:,n) = LF1(:,:,nn);
        else
            LF(:,:,n) = zeros(size(LF1(:,:,1)));
        end

    end
    clear LF1 L_reshaped
    if MFC % Account for mean field correction
        for n = 1:size(LF,3)
            LF(:,:,n) = S.M*LF(:,:,n);
        end
    end

    if do_SSS % Apply iterative SSS to lead fields
        ft_progress('init', 'text', 'Performing SSS correction...')      % ascii progress bar
        for n = 1:size(LF, 3)
            ft_progress(n/floor(size(LF,3)), 'Processing voxel %d of %d', n, floor(size(LF,3)))
            [sss_outs]=SSS_opm_lf(LF(:,:,n),sensor_info.pos,sensor_info.ors,S.Lin,S.Lout,S.sss_ins);
            LF(:,:,n) = sss_outs.phi_in;
        end
        ft_progress('close')
    end

    for n = size(sourcepos,1):-1:1
        % project to two tangential fields
        R = sourcepos(n,:) - origin;
        [etheta,ephi] = calctangent_all_bf(R);
        tan_ors = [etheta;ephi];
        ltan(:,:,n) = LF(:,:,n)*tan_ors';
    end

end

bf_outs.LF = LF;
bf_outs.ltan = ltan;
if do_bf
    ft_progress('init', 'text', 'Performing beamforming...')      % ascii progress bar
    for n = 1:size(sourcepos,1)
        % project to two tangential fields
        R = sourcepos(n,:) - origin;
        [etheta,ephi] = calctangent_all_bf(R);
        tan_ors = [etheta;ephi];
        ltan = LF(:,:,n)*tan_ors';

        % get optimal combination for max SNR
        [v,d] = svd(ltan'*S.Cinv*ltan);
        [~,id] = min(diag(d));
        lopt = ltan*v(:,id); % turn to nAm amplitude
        tan_or_opt = [tan_ors'*(v(:,id))]';

        bf_outs.lopt(:,:,n) = lopt;
        bf_outs.or_opt(n,:) = tan_or_opt;
        W1(:,n) = (S.Cinv*lopt/(lopt'*S.Cinv*lopt));
        Z{n} = [W1(:,n)'*S.C*W1(:,n)]./[W1(:,n)'*eye(size(W1(:,n),1))*W1(:,n)];
        if isfield(S,'Ca')
            Qa = W1(:,n)'*S.Ca*W1(:,n);
            Qc = W1(:,n)'*S.Cc*W1(:,n);
            T1{n} = (Qa-Qc)./(2*Qc);
        end
        ft_progress(n/size(sourcepos,1), 'Processing voxel %d of %d', n, size(sourcepos,1))
    end
    ft_progress('close')

    bf_outs.Weights = W1;
    bf_outs.Z_stat = Z;
    if isfield(S,'Ca')
        bf_outs.T_stat = T1;
        if strcmpi(max_or_min,'Max')
            cax_val = [0.5*max(cell2mat(T1)) max(cell2mat(T1))];
            [~, l] = max(cell2mat(T1));
        elseif strcmpi(max_or_min,'Min')
            cax_val = [min(cell2mat(T1)), 0.5.*min(cell2mat(T1))];
            [~, l] = min(cell2mat(T1));
        end
        try
            fg = figure;
            fg.Name = fwd_method;
            subplot(2,3,[1 2 4 5])
            ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none');
            hold on
            scatter3(sourcepos(:,1),sourcepos(:,2),sourcepos(:,3),50,cell2mat(T1),'filled');
            colormap hot;caxis([cax_val]);cb = colorbar;cb.Label.String = 'Tstat';
            axis equal
            view([120 20]);
            fig = gcf;
            fig.Color = 'w';
            ax = gca;
            ax.FontSize = 14;

            subplot(2,3,[3 6])
            % Find voxel with min/max value of Tstat
            dip_loc = sourcepos(l,:);
            dip_or = bf_outs.or_opt(l,:);
            ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none');
            hold on
            plot3(dip_loc(1),dip_loc(2),dip_loc(3),'.','markersize',20);
            quiver3(dip_loc(1),dip_loc(2),dip_loc(3),bf_outs.or_opt(l,1)./100,bf_outs.or_opt(l,2)./100,bf_outs.or_opt(l,3)./100);
            bf_outs.dip_loc = dip_loc;
            bf_outs.dip_or = dip_or;
        catch
            disp('No min or max T stat preference selected')
        end
    end
end

end

%%%%%%%%%%%% OTHER FUNCTIONS %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tanu, tanv] = calctangent_all_bf(RDip)
% Same as usual calctangent but put here for ease

x=RDip(1);
y=RDip(2);
z=RDip(3);
r=sqrt(x*x+y*y+z*z);

if (x==0) && (y==0)
    tanu(1)=1.0; tanu(2)=0; tanu(3)=0;
    tanv(1)=0; tanv(2)=1.0; tanv(3)=0;
else
    RZXY= -(r-z)*x*y;
    X2Y2= 1/(x*x+y*y);

    tanu(1)= (z*x*x + r*y*y) * X2Y2/r;
    tanu(2)= RZXY * X2Y2/r;
    tanu(3)= -x/r;

    tanv(1)= RZXY * X2Y2/r;
    tanv(2)= (z*y*y + r*x*x) * X2Y2/r;
    tanv(3)= -y/r;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Leads, L_reshaped] = Forward_OPM_BEM_stenroos_inscript(meshes_file,sensor_info,sourcepos)
% Use Mattis Stenroos code to create a BEM lead field model
% mri_file = .mri file
% sensor_info = sensor positions (Nx3), orientations(Nx3), and labels (Nx1)
% sourcepos = source positions

% Load OPM sensor positions
grad.coilpos = sensor_info.pos.*1;
grad.coilori = sensor_info.ors;
grad.label = sensor_info.label;
grad.chanpos = grad.coilpos;
grad.chanori = grad.coilori;
grad.units = 'm';

% Create head model
load(meshes_file)
% openmeeg_path = 'C:\Program Files\OpenMEEG\bin';
% if ~contains(path,openmeeg_path)
%     addpath(openmeeg_path)
% end
% curr_path = pwd;
% cd(openmeeg_path)

% Can get mesh intersection in the BEM model so redo the meshes to a
% smaller number of points (gets rid of spikes in meshes that can intersect)
cfg             = [];
cfg.tissue      = {'brain', 'skull', 'scalp'};
cfg.numvertices = [2500, 2500, 2500];
mesh = ft_prepare_mesh(cfg,segmentedmri);

% Add BEM framework paths. Put here your hbf root directory
tmp = matlab.desktop.editor.getActive;
func_code_dir = fileparts(tmp.Filename); % where the script is
hbfroot = [func_code_dir '\hbf_lc_p-master'];
thisdir = cd;
cd(hbfroot);
hbf_SetPaths();
cd(thisdir);

% 1. Boundary meshes
%
% Boundary meshes are described as a Nx1,cell array  where each cell is a
% struct that contains fields "p" for points (vertices), and "e" for
% element description (faces). The element description is 1-based.
% bmeshes{I} =
%
%   p: [Number of vertices x 3]
%   e: [Number of triangles x 3]

bmeshes = {};
bmeshes{1,1}.p = mesh(1).pos./1000;
bmeshes{1,1}.e = mesh(1).tri;
bmeshes{2,1}.p = mesh(2).pos./1000;
bmeshes{2,1}.e = mesh(2).tri;
bmeshes{3,1}.p = mesh(3).pos./1000;
bmeshes{3,1}.e = mesh(3).tri;

% First time with new meshing: tests...
% Problems with the BEM are typically due to unsuitable or ill-specified
% meshes. Even if meshes are otherwise good, the orientation of the
% triangles is often wrong; this is the most common reason for computations
% going wrong. This BEM framework assumes CCW orientation. To check and
% flip the orientation, you can just type

Nmeshes = length(bmeshes);
success = zeros(Nmeshes,1);
for M = 1:Nmeshes
    [bmeshes{M},success(M)] = hbf_CorrectTriangleOrientation(bmeshes{M});
end

% You can do further checks using, e.g., this simple tool. If orientation
% correction failed, this may give you at least some tips, where to look
% for the error.

status = cell(Nmeshes,1);
for M = 1:Nmeshes
    status{M} = hbf_CheckMesh(bmeshes{M});
end

% The convention of this BEM framework is that the innermost mesh has number
% 1 and outermost mesh M; this is mandatory. If the model is nested (like a
% 3-shell model), the order can be checked/set with
bmeshes = hbf_SortNestedMeshes(bmeshes);

figure
trisurf(bmeshes{1}.e,bmeshes{1}.p(:,1),bmeshes{1}.p(:,2),bmeshes{1}.p(:,3),...
    'FaceColor','r','EdgeColor','none','FaceAlpha',1)
hold on
trisurf(bmeshes{2}.e,bmeshes{2}.p(:,1),bmeshes{2}.p(:,2),bmeshes{2}.p(:,3),...
    'FaceColor','g','EdgeColor','none','FaceAlpha',0.5)
trisurf(bmeshes{3}.e,bmeshes{3}.p(:,1),bmeshes{3}.p(:,2),bmeshes{3}.p(:,3),...
    'FaceColor','b','EdgeColor','none','FaceAlpha',0.1)
axis equal
drawnow

% 2. MEG sensors
%
% Set MEG coils. There are some options for describing coils:
% coils =
%   QP: field computations points [Number of field computation points x 3]
%   QN: sensor orientation        [Number of field computation points x 3];
%       each must be unit length
%   QPinds: start and end indices of the QPs of each sensor [number of sensors x 2]
%   QW: sensor integration weights [number of field computation points x 1]
%   QtoC: integral weights and conversion QP -> sensors
%      [Number of sensors x Number of field computation points]
%   p:  sensor positions [Number of sensors x 3]
%   n:  sensor normals [Number of sensors x 3]; each must be unit-length.
%
% You must give (QP, QN and QtoC) or (QP, QN, QP and QPinds) or (p and n).
% If overlapping or conflicting information is given, the order of
% preference is the same as in this list.
%
% This sample coil set represents a 306-channel Neuromag system and uses
% 4-point sensor integrals as in, for example, Neuromag and MNE software. The
% coregistration is arbitrary, not the one in the data set. There is some
% further information in fields 'names','mesh3d', and 'description'.

coils = [];
coils.p = grad.chanpos;
coils.n = grad.chanori;

% 3. Source space
% We need to provide dipole locations. If we want to use oriented sources,
% we need to provide also the orientations.
%
% Here we can have a cortical mesh with also normal vectors for mesh
% vertices.
% cortex =
%
%   struct with fields:
%
%      p: [20484×3 single]
%      e: [40960×3 int32]
%     nn: [20484×3 double]
%
% so we can directly use cortex.p as source positions and cortex.nn as
% source orientations.

cortex.p = sourcepos;

% 4. Conductivities
% Set conductivities for the head model. In this example, we use a (nested)
% 3-layer model, so we need to give conductivities for the brain, skull,
% and scalp compartments. Remember to give the conductivities in the
% correct order!

ci=[1 1/50 1]*.33; % conductivity inside each surface --- remember the order!
co=[ci(2:3) 0];    % conductivity outside each surface


% ------------------------------------------------------------------------
% Now we should be ready for computation, unless something went wrong with
% the checks...

tic
% 1. BEM double-layer operators for potentials: do this once per meshing...
D = hbf_BEMOperatorsPhi_LC(bmeshes);

% 2. BEM operators for magnetic field due to volume currents:
% do this, if the boundary meshing or MEG sensors change

DB = hbf_BEMOperatorsB_Linear(bmeshes,coils);

% 3. Full transfer matrix
% Build BEM transfer matrix for potential using the isolated source
% approach, isolation set to surface 1.

Tphi_full = hbf_TM_Phi_LC_ISA2(D,ci,co,1);

% 4. Transfer matrices for volume component of the magnetic field and the
%potential

TBvol = hbf_TM_Bvol_Linear(DB,Tphi_full,ci,co);

% 5. Forward solutions = leadfield matrices
% LFM for orthogonal unit dipoles (orientations x,y,z in world coordinates).
% MEG only
Leads = hbf_LFM_B_LC(bmeshes,coils,TBvol,cortex.p);
fprintf('Total time was %ds.\n',round(toc));
L_reshaped = reshape(Leads,size(sensor_info.pos,1),3,[]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Lead_fields, L_reshaped] = Forward_OPM_BEM_FT_inscript(meshes_file,sensor_info,sourcepos)
% Use FieldTrip to create a BEM lead field model
% mri_file = .mri file
% sensor_info = sensor positions (Nx3), orientations(Nx3), and labels (Nx1)
% sourcepos = source positions

% Load OPM sensor positions
grad.coilpos = sensor_info.pos;
grad.coilori = sensor_info.ors;
grad.label = sensor_info.label;
grad.chanpos = grad.coilpos;
grad.chanori = grad.coilori;
grad.units = 'm';

% Create head model
load(meshes_file)
openmeeg_path = 'C:\Program Files\OpenMEEG\bin';
if ~contains(path,openmeeg_path)
    addpath(openmeeg_path)
end
curr_path = pwd;
cd(openmeeg_path)

% Dilate brain surface
% se = strel('sphere',1);
% segmentedmri.brain = imdilate(segmentedmri.brain,se);

% Can get mesh intersection in the BEM model so redo the meshes to a
% smaller number of points (gets rid of spikes in meshes that can intersect)
cfg             = [];
cfg.tissue      = {'brain', 'skull', 'scalp'};
cfg.numvertices = [2500, 500, 500];
mesh = ft_prepare_mesh(cfg,segmentedmri);
mesh = ft_convert_units(mesh,'m');

% figure
% ft_plot_mesh(mesh(1), 'facecolor','r', 'facealpha', 1, 'edgecolor', 'k', 'edgealpha', 1);
% hold on
% ft_plot_mesh(mesh(2), 'facecolor','g', 'facealpha', 0.4, 'edgecolor', 'k', 'edgealpha', 0.1);
% hold on
% ft_plot_mesh(mesh(3), 'facecolor','b', 'facealpha', 0.1, 'edgecolor', 'k', 'edgealpha', 0.1);

cfg = [];
cfg.method    = 'openmeeg';
vol = ft_prepare_headmodel(cfg, mesh);

% Visualise
headmodel = ft_convert_units(vol, 'm');

% Create leadfield
cfg                = [];
cfg.grad           = grad;
cfg.headmodel      = headmodel;
cfg.unit           = 'm';
cfg.sourcemodel.pos = sourcepos;
cfg.sourcemodel.inside = ones(size(sourcepos(:,1)));
% Remove points on mesh boundaries
% in = intriangulation_inscript(mesh(1).pos,mesh(1).tri,sourcepos);
% cfg.sourcemodel.inside(in) = 0;
% k = boundary(sourcepos,1);
% k = unique(k(:));
% cfg.sourcemodel.inside(k) = 0;

figure
ft_plot_sens(grad, 'style', '*b');
hold on
ft_plot_headmodel(headmodel, 'facecolor', 'cortex');
plot3(sourcepos(find(cfg.sourcemodel.inside),1),sourcepos(find(cfg.sourcemodel.inside),2),...
    sourcepos(find(cfg.sourcemodel.inside),3),'rx')

cfg.reducerank     = 2;
Lead_fields        = ft_prepare_leadfield(cfg);

inside_idx = find(Lead_fields.inside);
Leads = [];
disp('Extracting Lead Fields')
for n = 1:length(inside_idx)
    disp(100.*n./length(inside_idx))
    Leads = [Leads, Lead_fields.leadfield{inside_idx(n)}];
end
L_reshaped = reshape(Leads,size(sensor_info.pos,1),3,[]);
cd(curr_path)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Lead_fields, L_reshaped] = Forward_OPM_shell_FT_inscript(meshes_file,sensor_info,sourcepos)
% Use FieldTrip to create a single shell lead field model
% mri_file = .mri file
% sensor_info = sensor positions (Nx3), orientations(Nx3), and labels (Nx1)
% sourcepos = source positions

% Load OPM sensor positions
grad.coilpos = sensor_info.pos.*1;
grad.coilori = sensor_info.ors;
grad.label = sensor_info.label;
grad.chanpos = grad.coilpos;
grad.chanori = grad.coilori;
grad.units = 'm';

% Create head model
load(meshes_file)
cfg = [];
cfg.grad      = grad;
cfg.method    = 'singleshell';
cfg.tissue    = 'brain'; % will be constructed on the fly from white+grey+csf
vol = ft_prepare_headmodel(cfg, segmentedmri);

% Visualise
headmodel = ft_convert_units(vol, 'm');

figure
ft_plot_sens(grad, 'style', '*b');
hold on
ft_plot_headmodel(headmodel, 'facecolor', 'cortex');
plot3(sourcepos(:,1),sourcepos(:,2),sourcepos(:,3),'rx')

% Create leadfield
cfg                = [];
cfg.grad           = grad;
cfg.headmodel      = headmodel;
% cfg.resolution     = 0.4;
cfg.unit           = 'm';
% cfg.grid.pos       = sourcepos.*1;
cfg.sourcemodel.pos = sourcepos;
cfg.sourcemodel.inside = ones(size(sourcepos(:,1)));
cfg.reducerank     = 2;
Lead_fields        = ft_prepare_leadfield(cfg);

inside_idx = find(Lead_fields.inside);
Leads = [];
for n = 1:length(inside_idx)
    Leads = [Leads, Lead_fields.leadfield{inside_idx(n)}];
end
L_reshaped = reshape(Leads,size(sensor_info.pos,1),3,[]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Lead_fields, L_reshaped] = Forward_OPM_local_FT_inscript(meshes_file,sensor_info,sourcepos)
% Use FieldTrip to create a local spheres lead field model
% mri_file = .mri file
% sensor_info = sensor positions (Nx3), orientations(Nx3), and labels (Nx1)
% sourcepos = source positions

% Load OPM sensor positions
grad.coilpos = sensor_info.pos.*1;
grad.coilori = sensor_info.ors;
grad.label = sensor_info.label;
grad.chanpos = grad.coilpos;
grad.chanori = grad.coilori;
grad.units = 'm';

% Create head model
figure
load(meshes_file)
cfg = [];
cfg.grad      = grad;
cfg.method    = 'localspheres';
cfg.tissue    = 'brain'; % will be constructed on the fly from white+grey+csf
vol = ft_prepare_headmodel(cfg, segmentedmri);

% Visualise
headmodel = ft_convert_units(vol, 'm');

figure
ft_plot_sens(grad, 'style', '*b');
hold on
ft_plot_headmodel(headmodel, 'facecolor', 'cortex');
plot3(sourcepos(:,1),sourcepos(:,2),sourcepos(:,3),'rx')

% Create leadfield
cfg                = [];
cfg.grad           = grad;
cfg.headmodel      = headmodel;
% cfg.resolution     = 0.4;
cfg.unit           = 'm';
% cfg.grid.pos       = sourcepos.*1;
cfg.sourcemodel.pos = sourcepos;
cfg.sourcemodel.inside = ones(size(sourcepos(:,1)));
cfg.reducerank     = 2;
Lead_fields        = ft_prepare_leadfield(cfg);

inside_idx = find(Lead_fields.inside);
Leads = [];
for n = 1:length(inside_idx)
    Leads = [Leads, Lead_fields.leadfield{inside_idx(n)}];
end
L_reshaped = reshape(Leads,size(sensor_info.pos,1),3,[]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Lead_fields, L_reshaped] = Forward_OPM_single_FT_inscript(meshes_file,sensor_info,sourcepos)
% Use FieldTrip to create a single sphere lead field model
% mri_file = .mri file
% sensor_info = sensor positions (Nx3), orientations(Nx3), and labels (Nx1)
% sourcepos = source positions

% Load OPM sensor positions
grad.coilpos = sensor_info.pos.*1;
grad.coilori = sensor_info.ors;
grad.label = sensor_info.label;
grad.chanpos = grad.coilpos;
grad.chanori = grad.coilori;
grad.units = 'm';

% Create head model
load(meshes_file)
cfg = [];
cfg.grad      = grad;
cfg.method    = 'singlesphere';
cfg.tissue    = 'brain'; % will be constructed on the fly from white+grey+csf
vol = ft_prepare_headmodel(cfg, segmentedmri);

% Visualise
headmodel = ft_convert_units(vol, 'm');

figure
ft_plot_sens(grad, 'style', '*b');
hold on
ft_plot_headmodel(headmodel, 'facecolor', 'cortex');
plot3(sourcepos(:,1),sourcepos(:,2),sourcepos(:,3),'rx')

% Create leadfield
cfg                = [];
cfg.grad           = grad;
cfg.headmodel      = headmodel;
% cfg.resolution     = 0.4;
cfg.unit           = 'm';
% cfg.grid.pos       = sourcepos.*1;
cfg.sourcemodel.pos = sourcepos;
cfg.sourcemodel.inside = ones(size(sourcepos(:,1)));
cfg.reducerank     = 2;
Lead_fields        = ft_prepare_leadfield(cfg);

inside_idx = find(Lead_fields.inside);
Leads = [];
for n = 1:length(inside_idx)
    Leads = [Leads, Lead_fields.leadfield{inside_idx(n)}];
end
L_reshaped = reshape(Leads,size(sensor_info.pos,1),3,[]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function in = intriangulation_inscript(vertices,faces,testp,heavytest)
% intriangulation: Test points in 3d wether inside or outside a (closed) triangulation
% usage: in = intriangulation(vertices,faces,testp,heavytest)
%
% arguments: (input)
%  vertices   - points in 3d as matrix with three columns
%
%  faces      - description of triangles as matrix with three columns.
%               Each row contains three indices into the matrix of vertices
%               which gives the three cornerpoints of the triangle.
%
%  testp      - points in 3d as matrix with three columns
%
%  heavytest  - int n >= 0. Perform n additional randomized rotation tests.
%
% IMPORTANT: the set of vertices and faces has to form a watertight surface!
%
% arguments: (output)
%  in - a vector of length size(testp,1), containing 0 and 1.
%       in(nr) =  0: testp(nr,:) is outside the triangulation
%       in(nr) =  1: testp(nr,:) is inside the triangulation
%       in(nr) = -1: unable to decide for testp(nr,:)
%
% Thanks to Adam A for providing the FEX submission voxelise. The
% algorithms of voxelise form the algorithmic kernel of intriangulation.
%
% Thanks to Sven to discussions about speed and avoiding problems in
% special cases.
%
% Example usage:
%
%      n = 10;
%      vertices = rand(n, 3)-0.5; % Generate random points
%      tetra = delaunayn(vertices); % Generate delaunay triangulization
%      faces = freeBoundary(TriRep(tetra,vertices)); % use free boundary as triangulation
%      n = 1000;
%      testp = 2*rand(n,3)-1; % Generate random testpoints
%      in = intriangulation(vertices,faces,testp);
%      % Plot results
%      h = trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3));
%      set(h,'FaceColor','black','FaceAlpha',1/3,'EdgeColor','none');
%      hold on;
%      plot3(testp(:,1),testp(:,2),testp(:,3),'b.');
%      plot3(testp(in==1,1),testp(in==1,2),testp(in==1,3),'ro');
%
% See also: intetrahedron, tsearchn, inpolygon
%
% Author: Johannes Korsawe, heavily based on voxelise from Adam A.
% E-mail: johannes.korsawe@volkswagen.de
% Release: 1.3
% Release date: 25/09/2013
% check number of inputs
if nargin<3,
    fprintf('??? Error using ==> intriangulation\nThree input matrices are needed.\n');in=[];return;
end
if nargin==3,
    heavytest = 0;
end
% check size of inputs
if size(vertices,2)~=3 || size(faces,2)~=3 || size(testp,2)~=3,
    fprintf('??? Error using ==> intriagulation\nAll input matrices must have three columns.\n');in=[];return;
end
ipmax = max(faces(:));zerofound = ~isempty(find(faces(:)==0, 1));
if ipmax>size(vertices,1) || zerofound,
    fprintf('??? Error using ==> intriangulation\nThe triangulation data is defect. use trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3)) for test of deficiency.\n');return;
end
% loop for heavytest
inreturn = zeros(size(testp,1),1);VER = vertices;TESTP = testp;
for n = 1:heavytest+1,
    % Randomize
    if n>1,
        v=rand(1,3);D=rotmatrix(v/norm(v),rand*180/pi);vertices=VER*D;testp = TESTP*D;
    else,
        vertices=VER;
    end

    % Preprocessing data
    meshXYZ = zeros(size(faces,1),3,3);
    for loop = 1:3,
        meshXYZ(:,:,loop) = vertices(faces(:,loop),:);
    end
    % Basic idea (ingenious from FeX-submission voxelise):
    % If point is inside, it will cross the triangulation an uneven number of times in each direction (x, -x, y, -y, z, -z).

    % The function VOXELISEinternal is about 98% identical to its version inside voxelise.m.
    % This includes the elaborate comments. Thanks to Adam A!

    % z-direction:
    % intialization of results and correction list
    [in,cl] = VOXELISEinternal(testp(:,1),testp(:,2),testp(:,3),meshXYZ);

    % x-direction:
    % has only to be done for those points, that were not determinable in the first step --> cl
    [in2,cl2] = VOXELISEinternal(testp(cl,2),testp(cl,3),testp(cl,1),meshXYZ(:,[2,3,1],:));
    % Use results of x-direction that determined "inside"
    in(cl(in2==1)) = 1;
    % remaining indices with unclear result
    cl = cl(cl2);

    % y-direction:
    % has only to be done for those points, that were not determinable in the first and second step --> cl
    [in3,cl3] = VOXELISEinternal(testp(cl,3),testp(cl,1),testp(cl,2),meshXYZ(:,[3,1,2],:));

    % Use results of y-direction that determined "inside"
    in(cl(in3==1)) = 1;
    % remaining indices with unclear result
    cl = cl(cl3);

    % mark those indices, where all three tests have failed
    in(cl) = -1;

    if n==1,
        inreturn = in;  % Starting guess
    else,
        % if ALWAYS inside, use as inside!
        %        I = find(inreturn ~= in);
        %        inreturn(I(in(I)==0)) = 0;

        % if AT LEAST ONCE inside, use as inside!
        I = find(inreturn ~= in);
        inreturn(I(in(I)==1)) = 1;

    end
end
in = inreturn;
end
%==========================================================================
function [OUTPUT,correctionLIST] = VOXELISEinternal(testx,testy,testz,meshXYZ)
% Prepare logical array to hold the logical data:
OUTPUT = false(size(testx,1),1);
%Identify the min and max x,y coordinates of the mesh:
meshZmin = min(min(meshXYZ(:,3,:)));meshZmax = max(max(meshXYZ(:,3,:)));
%Identify the min and max x,y,z coordinates of each facet:
meshXYZmin = min(meshXYZ,[],3);meshXYZmax = max(meshXYZ,[],3);
%======================================================
% TURN OFF DIVIDE-BY-ZERO WARNINGS
%======================================================
%This prevents the Y1predicted, Y2predicted, Y3predicted and YRpredicted
%calculations creating divide-by-zero warnings.  Suppressing these warnings
%doesn't affect the code, because only the sign of the result is important.
%That is, 'Inf' and '-Inf' results are ok.
%The warning will be returned to its original state at the end of the code.
warningrestorestate = warning('query', 'MATLAB:divideByZero');
%warning off MATLAB:divideByZero
%======================================================
% START COMPUTATION
%======================================================
correctionLIST = [];   %Prepare to record all rays that fail the voxelisation.  This array is built on-the-fly, but since
%it ought to be relatively small should not incur too much of a speed penalty.
% Loop through each testpoint.
% The testpoint-array will be tested by passing rays in the z-direction through
% each x,y coordinate of the testpoints, and finding the locations where the rays cross the mesh.
facetCROSSLIST = zeros(1,1e3);  % uses countindex: nf
nm = size(meshXYZmin,1);
for loop = 1:length(OUTPUT),

    nf = 0;
    %    % - 1a - Find which mesh facets could possibly be crossed by the ray:
    %    possibleCROSSLISTy = find( meshXYZmin(:,2)<=testy(loop) & meshXYZmax(:,2)>=testy(loop) );
    %    % - 1b - Find which mesh facets could possibly be crossed by the ray:
    %    possibleCROSSLIST = possibleCROSSLISTy( meshXYZmin(possibleCROSSLISTy,1)<=testx(loop) & meshXYZmax(possibleCROSSLISTy,1)>=testx(loop) );
    % Do - 1a - and - 1b - faster
    possibleCROSSLISTy = find((testy(loop)-meshXYZmin(:,2)).*(meshXYZmax(:,2)-testy(loop))>0);
    possibleCROSSLISTx = (testx(loop)-meshXYZmin(possibleCROSSLISTy,1)).*(meshXYZmax(possibleCROSSLISTy,1)-testx(loop))>0;
    possibleCROSSLIST = possibleCROSSLISTy(possibleCROSSLISTx);

    if isempty(possibleCROSSLIST)==0  %Only continue the analysis if some nearby facets were actually identified

        % - 2 - For each facet, check if the ray really does cross the facet rather than just passing it close-by:

        % GENERAL METHOD:
        % 1. Take each edge of the facet in turn.
        % 2. Find the position of the opposing vertex to that edge.
        % 3. Find the position of the ray relative to that edge.
        % 4. Check if ray is on the same side of the edge as the opposing vertex.
        % 5. If this is true for all three edges, then the ray definitely passes through the facet.
        %
        % NOTES:
        % 1. If the ray crosses exactly on an edge, this is counted as crossing the facet.
        % 2. If a ray crosses exactly on a vertex, this is also taken into account.

        for loopCHECKFACET = possibleCROSSLIST'

            %Check if ray crosses the facet.  This method is much (>>10 times) faster than using the built-in function 'inpolygon'.
            %Taking each edge of the facet in turn, check if the ray is on the same side as the opposing vertex.  If so, let testVn=1

            Y1predicted = meshXYZ(loopCHECKFACET,2,2) - ((meshXYZ(loopCHECKFACET,2,2)-meshXYZ(loopCHECKFACET,2,3)) * (meshXYZ(loopCHECKFACET,1,2)-meshXYZ(loopCHECKFACET,1,1))/(meshXYZ(loopCHECKFACET,1,2)-meshXYZ(loopCHECKFACET,1,3)));
            YRpredicted = meshXYZ(loopCHECKFACET,2,2) - ((meshXYZ(loopCHECKFACET,2,2)-meshXYZ(loopCHECKFACET,2,3)) * (meshXYZ(loopCHECKFACET,1,2)-testx(loop))/(meshXYZ(loopCHECKFACET,1,2)-meshXYZ(loopCHECKFACET,1,3)));

            if (Y1predicted > meshXYZ(loopCHECKFACET,2,1) && YRpredicted > testy(loop)) || (Y1predicted < meshXYZ(loopCHECKFACET,2,1) && YRpredicted < testy(loop)) || (meshXYZ(loopCHECKFACET,2,2)-meshXYZ(loopCHECKFACET,2,3)) * (meshXYZ(loopCHECKFACET,1,2)-testx(loop)) == 0
                %                testV1 = 1;   %The ray is on the same side of the 2-3 edge as the 1st vertex.
            else
                %                testV1 = 0;   %The ray is on the opposite side of the 2-3 edge to the 1st vertex.
                % As the check is for ALL three checks to be true, we can continue here, if only one check fails
                continue;
            end %if

            Y2predicted = meshXYZ(loopCHECKFACET,2,3) - ((meshXYZ(loopCHECKFACET,2,3)-meshXYZ(loopCHECKFACET,2,1)) * (meshXYZ(loopCHECKFACET,1,3)-meshXYZ(loopCHECKFACET,1,2))/(meshXYZ(loopCHECKFACET,1,3)-meshXYZ(loopCHECKFACET,1,1)));
            YRpredicted = meshXYZ(loopCHECKFACET,2,3) - ((meshXYZ(loopCHECKFACET,2,3)-meshXYZ(loopCHECKFACET,2,1)) * (meshXYZ(loopCHECKFACET,1,3)-testx(loop))/(meshXYZ(loopCHECKFACET,1,3)-meshXYZ(loopCHECKFACET,1,1)));
            if (Y2predicted > meshXYZ(loopCHECKFACET,2,2) && YRpredicted > testy(loop)) || (Y2predicted < meshXYZ(loopCHECKFACET,2,2) && YRpredicted < testy(loop)) || (meshXYZ(loopCHECKFACET,2,3)-meshXYZ(loopCHECKFACET,2,1)) * (meshXYZ(loopCHECKFACET,1,3)-testx(loop)) == 0
                %                testV2 = 1;   %The ray is on the same side of the 3-1 edge as the 2nd vertex.
            else
                %                testV2 = 0;   %The ray is on the opposite side of the 3-1 edge to the 2nd vertex.
                % As the check is for ALL three checks to be true, we can continue here, if only one check fails
                continue;
            end %if

            Y3predicted = meshXYZ(loopCHECKFACET,2,1) - ((meshXYZ(loopCHECKFACET,2,1)-meshXYZ(loopCHECKFACET,2,2)) * (meshXYZ(loopCHECKFACET,1,1)-meshXYZ(loopCHECKFACET,1,3))/(meshXYZ(loopCHECKFACET,1,1)-meshXYZ(loopCHECKFACET,1,2)));
            YRpredicted = meshXYZ(loopCHECKFACET,2,1) - ((meshXYZ(loopCHECKFACET,2,1)-meshXYZ(loopCHECKFACET,2,2)) * (meshXYZ(loopCHECKFACET,1,1)-testx(loop))/(meshXYZ(loopCHECKFACET,1,1)-meshXYZ(loopCHECKFACET,1,2)));
            if (Y3predicted > meshXYZ(loopCHECKFACET,2,3) && YRpredicted > testy(loop)) || (Y3predicted < meshXYZ(loopCHECKFACET,2,3) && YRpredicted < testy(loop)) || (meshXYZ(loopCHECKFACET,2,1)-meshXYZ(loopCHECKFACET,2,2)) * (meshXYZ(loopCHECKFACET,1,1)-testx(loop)) == 0
                %                testV3 = 1;   %The ray is on the same side of the 1-2 edge as the 3rd vertex.
            else
                %                testV3 = 0;   %The ray is on the opposite side of the 1-2 edge to the 3rd vertex.
                % As the check is for ALL three checks to be true, we can continue here, if only one check fails
                continue;
            end %if

            nf=nf+1;facetCROSSLIST(nf)=loopCHECKFACET;

        end %for
        % Use only values ~=0
        facetCROSSLIST = facetCROSSLIST(1:nf);

        % - 3 - Find the z coordinate of the locations where the ray crosses each facet:
        gridCOzCROSS = zeros(1,nf);
        for loopFINDZ = facetCROSSLIST

            % METHOD:
            % 1. Define the equation describing the plane of the facet.  For a
            % more detailed outline of the maths, see:
            % http://local.wasp.uwa.edu.au/~pbourke/geometry/planeeq/
            %    Ax + By + Cz + D = 0
            %    where  A = y1 (z2 - z3) + y2 (z3 - z1) + y3 (z1 - z2)
            %           B = z1 (x2 - x3) + z2 (x3 - x1) + z3 (x1 - x2)
            %           C = x1 (y2 - y3) + x2 (y3 - y1) + x3 (y1 - y2)
            %           D = - x1 (y2 z3 - y3 z2) - x2 (y3 z1 - y1 z3) - x3 (y1 z2 - y2 z1)
            % 2. For the x and y coordinates of the ray, solve these equations to find the z coordinate in this plane.

            planecoA = meshXYZ(loopFINDZ,2,1)*(meshXYZ(loopFINDZ,3,2)-meshXYZ(loopFINDZ,3,3)) + meshXYZ(loopFINDZ,2,2)*(meshXYZ(loopFINDZ,3,3)-meshXYZ(loopFINDZ,3,1)) + meshXYZ(loopFINDZ,2,3)*(meshXYZ(loopFINDZ,3,1)-meshXYZ(loopFINDZ,3,2));
            planecoB = meshXYZ(loopFINDZ,3,1)*(meshXYZ(loopFINDZ,1,2)-meshXYZ(loopFINDZ,1,3)) + meshXYZ(loopFINDZ,3,2)*(meshXYZ(loopFINDZ,1,3)-meshXYZ(loopFINDZ,1,1)) + meshXYZ(loopFINDZ,3,3)*(meshXYZ(loopFINDZ,1,1)-meshXYZ(loopFINDZ,1,2));
            planecoC = meshXYZ(loopFINDZ,1,1)*(meshXYZ(loopFINDZ,2,2)-meshXYZ(loopFINDZ,2,3)) + meshXYZ(loopFINDZ,1,2)*(meshXYZ(loopFINDZ,2,3)-meshXYZ(loopFINDZ,2,1)) + meshXYZ(loopFINDZ,1,3)*(meshXYZ(loopFINDZ,2,1)-meshXYZ(loopFINDZ,2,2));
            planecoD = - meshXYZ(loopFINDZ,1,1)*(meshXYZ(loopFINDZ,2,2)*meshXYZ(loopFINDZ,3,3)-meshXYZ(loopFINDZ,2,3)*meshXYZ(loopFINDZ,3,2)) - meshXYZ(loopFINDZ,1,2)*(meshXYZ(loopFINDZ,2,3)*meshXYZ(loopFINDZ,3,1)-meshXYZ(loopFINDZ,2,1)*meshXYZ(loopFINDZ,3,3)) - meshXYZ(loopFINDZ,1,3)*(meshXYZ(loopFINDZ,2,1)*meshXYZ(loopFINDZ,3,2)-meshXYZ(loopFINDZ,2,2)*meshXYZ(loopFINDZ,3,1));

            if abs(planecoC) < 1e-14
                planecoC=0;
            end

            gridCOzCROSS(facetCROSSLIST==loopFINDZ) = (- planecoD - planecoA*testx(loop) - planecoB*testy(loop)) / planecoC;

        end %for
        if isempty(gridCOzCROSS),continue;end

        %Remove values of gridCOzCROSS which are outside of the mesh limits (including a 1e-12 margin for error).
        gridCOzCROSS = gridCOzCROSS( gridCOzCROSS>=meshZmin-1e-12 & gridCOzCROSS<=meshZmax+1e-12 );
        %Round gridCOzCROSS to remove any rounding errors, and take only the unique values:
        gridCOzCROSS = round(gridCOzCROSS*1e10)/1e10;

        % Replacement of the call to unique (gridCOzCROSS = unique(gridCOzCROSS);) by the following line:
        tmp = sort(gridCOzCROSS);I=[0,tmp(2:end)-tmp(1:end-1)]~=0;gridCOzCROSS = [tmp(1),tmp(I)];

        % - 4 - Label as being inside the mesh all the voxels that the ray passes through after crossing one facet before crossing another facet:

        if rem(numel(gridCOzCROSS),2)==0  % Only rays which cross an even number of facets are voxelised

            for loopASSIGN = 1:(numel(gridCOzCROSS)/2)
                voxelsINSIDE = (testz(loop)>gridCOzCROSS(2*loopASSIGN-1) & testz(loop)<gridCOzCROSS(2*loopASSIGN));
                OUTPUT(loop) = voxelsINSIDE;
                if voxelsINSIDE,break;end
            end %for

        elseif numel(gridCOzCROSS)~=0    % Remaining rays which meet the mesh in some way are not voxelised, but are labelled for correction later.
            correctionLIST = [ correctionLIST; loop ];
        end %if

    end %if

end %for
%======================================================
% RESTORE DIVIDE-BY-ZERO WARNINGS TO THE ORIGINAL STATE
%======================================================
warning(warningrestorestate)
% J.Korsawe: A correction is not possible as the testpoints need not to be
%            ordered in any way.
%            voxelise contains a correction algorithm which is appended here
%            without changes in syntax.
return
%======================================================
% USE INTERPOLATION TO FILL IN THE RAYS WHICH COULD NOT BE VOXELISED
%======================================================
%For rays where the voxelisation did not give a clear result, the ray is
%computed by interpolating from the surrounding rays.
countCORRECTIONLIST = size(correctionLIST,1);
if countCORRECTIONLIST>0

    %If necessary, add a one-pixel border around the x and y edges of the
    %array.  This prevents an error if the code tries to interpolate a ray at
    %the edge of the x,y grid.
    if min(correctionLIST(:,1))==1 || max(correctionLIST(:,1))==numel(gridCOx) || min(correctionLIST(:,2))==1 || max(correctionLIST(:,2))==numel(gridCOy)
        gridOUTPUT     = [zeros(1,voxcountY+2,voxcountZ);zeros(voxcountX,1,voxcountZ),gridOUTPUT,zeros(voxcountX,1,voxcountZ);zeros(1,voxcountY+2,voxcountZ)];
        correctionLIST = correctionLIST + 1;
    end

    for loopC = 1:countCORRECTIONLIST
        voxelsforcorrection = squeeze( sum( [ gridOUTPUT(correctionLIST(loopC,1)-1,correctionLIST(loopC,2)-1,:) ,...
            gridOUTPUT(correctionLIST(loopC,1)-1,correctionLIST(loopC,2),:)   ,...
            gridOUTPUT(correctionLIST(loopC,1)-1,correctionLIST(loopC,2)+1,:) ,...
            gridOUTPUT(correctionLIST(loopC,1),correctionLIST(loopC,2)-1,:)   ,...
            gridOUTPUT(correctionLIST(loopC,1),correctionLIST(loopC,2)+1,:)   ,...
            gridOUTPUT(correctionLIST(loopC,1)+1,correctionLIST(loopC,2)-1,:) ,...
            gridOUTPUT(correctionLIST(loopC,1)+1,correctionLIST(loopC,2),:)   ,...
            gridOUTPUT(correctionLIST(loopC,1)+1,correctionLIST(loopC,2)+1,:) ,...
            ] ) );
        voxelsforcorrection = (voxelsforcorrection>=4);
        gridOUTPUT(correctionLIST(loopC,1),correctionLIST(loopC,2),voxelsforcorrection) = 1;
    end %for
    %Remove the one-pixel border surrounding the array, if this was added
    %previously.
    if size(gridOUTPUT,1)>numel(gridCOx) || size(gridOUTPUT,2)>numel(gridCOy)
        gridOUTPUT = gridOUTPUT(2:end-1,2:end-1,:);
    end

end %if
%disp([' Ray tracing result: ',num2str(countCORRECTIONLIST),' rays (',num2str(countCORRECTIONLIST/(voxcountX*voxcountY)*100,'%5.1f'),'% of all rays) exactly crossed a facet edge and had to be computed by interpolation.'])
end %function
%==========================================================================
function D = rotmatrix(v,deg)
% calculate the rotation matrix about v by deg degrees
deg=deg/180*pi;
if deg~=0,
    v=v/norm(v);
    v1=v(1);v2=v(2);v3=v(3);ca=cos(deg);sa=sin(deg);
    D=[ca+v1*v1*(1-ca),v1*v2*(1-ca)-v3*sa,v1*v3*(1-ca)+v2*sa;
        v2*v1*(1-ca)+v3*sa,ca+v2*v2*(1-ca),v2*v3*(1-ca)-v1*sa;
        v3*v1*(1-ca)-v2*sa,v3*v2*(1-ca)+v1*sa,ca+v3*v3*(1-ca)];
else,
    D=eye(3,3);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%