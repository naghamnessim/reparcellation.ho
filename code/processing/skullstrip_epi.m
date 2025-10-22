function [mean_epi_stripped_path, epi_stripped_path] = skullstrip_epi( ...
    realigned_epi_path, brain_mask_epi_path, t1_stripped_path, varargin)
% Inputs:
%   realigned_epi_path  (char) — motion-corrected 4D EPI (EPI→T1 header already updated)
%   brain_mask_epi_path (char) — brain mask already on the EPI grid (nearest-neighbor resampled)
%   t1_stripped_path    (char) — skull-stripped T1 (used only for QC overlay)
%   Name-Value:
%     'SaveFigures' (logical, default false) — save QC PNGs
%     'SaveDir'     (char, default '')       — folder to save QC PNGs (created if needed)
%     'QCOnly'      (logical, default false) — do not recompute data; only regenerate QC
%
% Outputs:
%   mean_epi_stripped_path (char) — path to skull-stripped mean EPI (3D)
%   epi_stripped_path      (char) — path to skull-stripped 4D EPI
%
% Purpose:
%   Apply an EPI-space brain mask (binary label image) to both the mean EPI and
%   the full 4D EPI, then generate light QC figures (pre/post/overlay). A single,
%   consistent mask is used voxelwise across all timepoints to avoid time-varying
%   censoring artifacts.
%
% Author: Nagham Nessim
% University of Geneva, 2025


p = inputParser;
addParameter(p,'SaveFigures',false,@islogical);
addParameter(p,'SaveDir','',@(s)ischar(s)||isstring(s));
addParameter(p,'QCOnly',false,@islogical);                                                        %  QC-only mode reads output and just rebuilds figures
parse(p,varargin{:});
opt = p.Results;

if opt.SaveFigures && ~isempty(opt.SaveDir) && ~exist(opt.SaveDir,'dir')
    mkdir(opt.SaveDir);
end

% -----------------------
% output paths
% -----------------------
outdir                 = fileparts(realigned_epi_path);
mean_epi_path          = fullfile(outdir, 'mean_epi.nii');
mean_epi_stripped_path = fullfile(outdir, 'mean_epi_stripped.nii');
epi_stripped_path      = fullfile(outdir, ['stripped_' spm_file(realigned_epi_path,'filename')]);
tr_path                = fullfile(outdir, ['tr' spm_file(mean_epi_stripped_path,'filename')]);


% -----------------------
% QC mode (no recomputation)
% -----------------------
if opt.QCOnly
    % Ensure mean_epi exists (quick compute if missing)
    if ~isfile(mean_epi_path)
        epi_vol  = spm_vol(realigned_epi_path);
        epi_data = spm_read_vols(epi_vol);
        mean_epi = mean(epi_data,4);
        mv = epi_vol(1); mv.fname = mean_epi_path;
        spm_write_vol(mv, mean_epi);
    else
        mean_epi = spm_read_vols(spm_vol(mean_epi_path));                                         % Reuse existing mean if present
    end

    % Ensure stripped mean exists
    if ~isfile(mean_epi_stripped_path)
        error('QCOnly requested but missing: %s', mean_epi_stripped_path);
    end

    % Ensure a mean-EPI to T1 reslice exists for overlay figure
    % Executes reslice for visualization only
    if ~isfile(tr_path)
        spm('defaults','FMRI'); spm_jobman('initcfg');
        matlabbatch = [];    
        matlabbatch{1}.spm.spatial.coreg.write.ref                 = { [t1_stripped_path ',1'] };            % set ref image as T1 stripped (target grid)
        matlabbatch{1}.spm.spatial.coreg.write.source              = { [mean_epi_stripped_path ',1'] };      % set source as mean stripped EPI 
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp     = 4;                                      % 4th-degree B-spline interpolation (continuous data)
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap       = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.write.roptions.mask       = 0;
        matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix     = 'tr';                                   % 'tr' prefix
        spm_jobman('run', matlabbatch);
    end

    % QC (Quality Check) figures (pre-stripping, post-stripping, overlay)
    try
        f1 = figure('Color','w'); imagesc(squeeze(mean_epi(:,:,round(end/2)))); axis image off
        title('Mean EPI (pre-strip)'); colormap(gray);
    if opt.SaveFigures && ~isempty(opt.SaveDir)
        png_out = fullfile(opt.SaveDir,'step5_mean_epi_pre_strip.png');
        fig_out = fullfile(opt.SaveDir,'step5_mean_epi_pre_strip.fig');
        exportgraphics(f1, png_out, 'Resolution', 200);
        savefig(f1, fig_out);
        fprintf('Saved (Pre-strip) in:\n  %s\n  %s\n', png_out, fig_out);
    end
    catch ME
        warning(ME.identifier, 'Step 5 pre-strip figure not saved: %s', ME.message);
    end

    try
        vol_post = spm_read_vols(spm_vol(mean_epi_stripped_path));
        f2 = figure('Color','w'); imagesc(squeeze(vol_post(:,:,round(end/2)))); axis image off
        title('Mean EPI (post-strip)'); colormap(gray);
        if opt.SaveFigures && ~isempty(opt.SaveDir)
            png_out = fullfile(opt.SaveDir,'step5_mean_epi_post_strip.png');
            fig_out = fullfile(opt.SaveDir,'step5_mean_epi_post_strip.fig');
            exportgraphics(f2, png_out, 'Resolution', 200);
            savefig(f2, fig_out);
            fprintf('Saved (Post-strip)in:\n  %s\n  %s\n', png_out, fig_out);
        end
    catch ME
        warning(ME.identifier, 'Step 5 post-strip figure not saved: %s', ME.message);
    end

    try
        t1_img   = spm_read_vols(spm_vol(t1_stripped_path));
        resl_epi = spm_read_vols(spm_vol(tr_path));
        k        = round(size(t1_img,3)/2);
        f3 = figure('Color','w');
        imshowpair(mat2gray(t1_img(:,:,k)), mat2gray(resl_epi(:,:,k)), 'falsecolor');
        title('T1 and Functional Image Overlay (Axial Slice)');
        if opt.SaveFigures && ~isempty(opt.SaveDir)
            png_out = fullfile(opt.SaveDir,'step5_overlay.png');
            fig_out = fullfile(opt.SaveDir,'step5_overlay.fig');
            exportgraphics(f3, png_out, 'Resolution', 200);
            savefig(f3, fig_out);
            fprintf('Saved Overlay in:\n  %s\n  %s\n', png_out, fig_out);
        end
    catch ME
        warning(ME.identifier, 'Step 5 overlay figure not saved: %s', ME.message);
    end

    return
end


% -----------------------
% Computation then visulaization
% -----------------------
% Mean EPI (3D) from the 4D
epi_vol  = spm_vol(realigned_epi_path);
epi_data = spm_read_vols(epi_vol);           % Read entire EPI array [X Y Z T]
mean_epi = mean(epi_data,4);                 % [X Y Z]
mv = epi_vol(1); mv.fname = mean_epi_path;
spm_write_vol(mv, mean_epi);

% Load the already-resliced EPI-space mask (binary)
mask = spm_read_vols(spm_vol(brain_mask_epi_path)) > 0.5;   % 3D binary mask (mask was resampled using Nearest neighbor interpolation to keep it binary)
mask = double(mask);

% Apply mask to mean EPI (3D)
mv.fname = mean_epi_stripped_path;                          % Output path for stripped mean image.
spm_write_vol(mv, mean_epi .* mask);                        % apply mask using element-wise multiplication

% Apply mask to full 4D EPI 
M = epi_vol(1).mat;                                         % Enforce consistent affine transformation across all frames
for i = 2:numel(epi_vol)
    epi_vol(i).mat = M;                                     % Normalize orientation of all volumes to the first volume
end
epi_stripped_4d = bsxfun(@times, epi_data, mask);           % Apply same 3D mask to all T timepoints
for t = 1:numel(epi_vol)
    epi_vol(t).fname = epi_stripped_path;                   % Write all frames into a single 4D EPI
    epi_vol(t).n     = [t 1];
    spm_write_vol(epi_vol(t), epi_stripped_4d(:,:,:,t));
end

% Reslice stripped mean EPI to T1 for overlay figure
spm('defaults','FMRI'); spm_jobman('initcfg');
matlabbatch = [];
matlabbatch{1}.spm.spatial.coreg.write.ref                 = { [t1_stripped_path ',1'] };
matlabbatch{1}.spm.spatial.coreg.write.source              = { [mean_epi_stripped_path ',1'] };
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp     = 4;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap       = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask       = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix     = 'tr';
spm_jobman('run', matlabbatch);

% QC figures (pre, post, overlay)
try
    f1 = figure('Color','w'); imagesc(squeeze(mean_epi(:,:,round(end/2)))); axis image off
    title('Mean EPI (pre-strip)'); colormap(gray);
    if opt.SaveFigures && ~isempty(opt.SaveDir)
        png_out = fullfile(opt.SaveDir,'step5_mean_epi_pre_strip.png');
        fig_out = fullfile(opt.SaveDir,'step5_mean_epi_pre_strip.fig');
        exportgraphics(f1, png_out, 'Resolution', 200);
        savefig(f1, fig_out);
        fprintf('Saved (Pre-strip) in:\n  %s\n  %s\n', png_out, fig_out);
    end

catch ME
    warning(ME.identifier, 'Step 5 pre-strip figure not saved: %s', ME.message);
end

try
    vol_post = spm_read_vols(spm_vol(mean_epi_stripped_path));
    f2 = figure('Color','w'); imagesc(squeeze(vol_post(:,:,round(end/2)))); axis image off
    title('Mean EPI (post-strip)'); colormap(gray);
    if opt.SaveFigures && ~isempty(opt.SaveDir)
        png_out = fullfile(opt.SaveDir,'step5_mean_epi_post_strip.png');
        fig_out = fullfile(opt.SaveDir,'step5_mean_epi_post_strip.fig');
        exportgraphics(f2, png_out, 'Resolution', 200);
        savefig(f2, fig_out);
        fprintf('Saved (Post-strip) in:\n  %s\n  %s\n', png_out, fig_out);
    end

catch ME
    warning(ME.identifier, 'Step 5 post-strip figure not saved: %s', ME.message);
end

try
    t1_img   = spm_read_vols(spm_vol(t1_stripped_path));
    resl_epi = spm_read_vols(spm_vol(tr_path));
    k        = round(size(t1_img,3)/2);
    f3 = figure('Color','w');
    imshowpair(mat2gray(t1_img(:,:,k)), mat2gray(resl_epi(:,:,k)), 'falsecolor');
    title('T1 and Functional Image Overlay (Axial Slice)');
    if opt.SaveFigures && ~isempty(opt.SaveDir)
            png_out = fullfile(opt.SaveDir,'step5_overlay.png');
            fig_out = fullfile(opt.SaveDir,'step5_overlay.fig');
            exportgraphics(f3, png_out, 'Resolution', 200);
            savefig(f3, fig_out);
            fprintf('Saved Overlay in:\n  %s\n  %s\n', png_out, fig_out);
    end
catch ME
    warning(ME.identifier, 'Step 5 overlay figure not saved: %s', ME.message);
end
end

