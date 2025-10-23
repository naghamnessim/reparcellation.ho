function [mean_t1_stripped, brain_mask_path] = skullstrip_T1(t1_path_1, t1_path_2, varargin)
% Inputs:
%   t1_path_1  (char) — path to first T1
%   t1_path_2  (char) — path to second T1
%   Name-Value:
%     'SaveFigures'  (logical, default false) — save QC PNGs
%     'SaveDir'      (char, default '')       — folder to save QC PNGs
%     'FigurePrefix' (char, default 'step2_') — prefix for QC PNG filenames
%     'QCOnly'       (logical, default false) — don’t recompute; only build Quality Check (QC) figs from existing outputs
%
% Outputs:
%   mean_t1_stripped (char) — path to skull-stripped mean T1
%   brain_mask_path  (char) — path to brain mask (GM+WM+CSF thresholded) in mean-T1 space
%
% Purpose:
%   Segment T1 images using SPM unified segmentation (bias field + tissue priors + MoG),
%   derive a robust intracranial mask from GM/WM/CSF, apply it to the (bias-corrected) mean
%   of two T1 runs, and generate QC figures.
%
% Author: Nagham Nessim
% University of Genova, 2025


p = inputParser;                                                                     
addParameter(p,'SaveFigures',false,@islogical);                                      % Toggle figure saving to disk.
addParameter(p,'SaveDir','',@(s)ischar(s)||isstring(s));                             % Output directory for QC.
addParameter(p,'FigurePrefix','step2_',@(s)ischar(s)||isstring(s));                  % Prefix for the result photo
addParameter(p,'QCOnly',false,@islogical);                                           % QC-only mode prevents recomputation which saves time.
parse(p,varargin{:});
opt = p.Results;

if opt.SaveFigures && ~isempty(opt.SaveDir) && ~exist(opt.SaveDir,'dir')
    mkdir(opt.SaveDir);
end
prefix = char(opt.FigurePrefix);
if ~isempty(prefix) && ~endsWith(prefix,'_'), prefix = [prefix '_']; end

% Keep all produced NIfTIs beside the first T1
outdir = fileparts(t1_path_1);

% Convenience paths that this function writes (or expects in QCOnly)
m1_path            = fullfile(outdir, ['m' spm_file(t1_path_1,'filename')]);
mean_t1_path       = fullfile(outdir, 'mean_T1.nii');
brain_mask_path    = fullfile(outdir, 'brain_mask_meanT1.nii');
mean_t1_stripped   = fullfile(outdir, 'mean_T1_stripped.nii');
single_stripped_fn = fullfile(outdir, 'T1_brain_stripped.nii');

% -----------------------
% QC mode only (no recomputation)
% -----------------------
if opt.QCOnly
    if ~isfile(mean_t1_stripped) || ~isfile(brain_mask_path)
        error('skullstrip_T1:QCOnlyMissing', ...
              'QCOnly requested but required files missing:\n%s\n%s', ...
              mean_t1_stripped, brain_mask_path);                                             % Fail fast if products aren’t present; QC-only should never recompute
    end

    %2x2 comparison panel
    try
        mean_orig        = spm_read_vols(spm_vol(mean_t1_path));
        orig_t1_biascorr = spm_read_vols(spm_vol(m1_path));
        single_stripped  = spm_read_vols(spm_vol(single_stripped_fn));
        stripped_data    = spm_read_vols(spm_vol(mean_t1_stripped));
        sl = round(size(orig_t1_biascorr,3)/2);                                               % Middle axial slice—quick, consistent QC position.

        f2 = figure('Name','Comparison of Original vs Skull-Stripped T1','Color','w');
        sgtitle('Comparison of Original vs Skull-Stripped T1');
        subplot(2,2,1); imshow(mean_orig(:,:,sl), []);        title('Mean T1');
        subplot(2,2,2); imshow(stripped_data(:,:,sl), []);    title('Mean T1 (stripped)');
        subplot(2,2,3); imshow(orig_t1_biascorr(:,:,sl), []); title('T1#1 (bias-corrected)');
        subplot(2,2,4); imshow(single_stripped(:,:,sl), []);  title('T1#1 (stripped)');

        if opt.SaveFigures && ~isempty(opt.SaveDir)
            png_out = fullfile(opt.SaveDir,[prefix 'T1_comparison.png']);
            fig_out = fullfile(opt.SaveDir,[prefix 'T1_comparison.fig']);
            exportgraphics(f2, png_out, 'Resolution', 200);
            savefig(f2, fig_out);
            fprintf('Saved:\n  %s\n  %s\n', png_out, fig_out);
        end
    catch ME
        warning(ME.identifier,'Step 2 comparison figure not saved: %s', ME.message);
    end
    return;  % IMPORTANT: do not recompute when QCOnly=true
end


% -----------------------
% Compute then visualize
% -----------------------
% Handle .gz on T1#2
if endsWith(t1_path_2, '.gz', 'IgnoreCase', true)
    files = gunzip(t1_path_2);                                                               % returns cell array of unzipped files
    t1_path_2 = files{1};                                                                    % Use unzipped path for subsequent operations.
end

% A) Segment T1#1  -> c1/c2/c3 + bias-corrected m*
spm('defaults','FMRI'); spm_jobman('initcfg');
matlabbatch = [];
matlabbatch{1}.spm.spatial.preproc.channel.vols     = {[t1_path_1 ',1']};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg  = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write    = [0 1];                                 % Save bias-corrected image (m*)

%SPM’s unified segmentation models intensity with mixtures of Gaussians and uses tissue priors (TPMs) plus bias-field correction to estimate GM/WM/CSF; 
%You only need GM+WM+CSF (classes 1–3) to form a reliable brain mask; 
%other classes (4–6) describe non-brain structures you intend to remove, so you don’t bother writing them.

for t = 1:3                                                       % Tissue priors GM/WM/CSF (prob maps).
    matlabbatch{1}.spm.spatial.preproc.tissue(t).tpm    = {fullfile(spm('Dir'),'tpm',['TPM.nii,' num2str(t)])};
    matlabbatch{1}.spm.spatial.preproc.tissue(t).ngaus  = 1;      % 1 Gaussian per tissue class
    matlabbatch{1}.spm.spatial.preproc.tissue(t).native = [1 0];  % write native space to build mask
    matlabbatch{1}.spm.spatial.preproc.tissue(t).warped = [0 0];  % no warping 
end

for t = 4:6                                                       % non-brain tissues
    matlabbatch{1}.spm.spatial.preproc.tissue(t).tpm    = {fullfile(spm('Dir'),'tpm',['TPM.nii,' num2str(t)])};
    matlabbatch{1}.spm.spatial.preproc.tissue(t).ngaus  = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(t).native = [0 0];  % do not write maps for these tissue types (won't be used for mask)
    matlabbatch{1}.spm.spatial.preproc.tissue(t).warped = [0 0];
end
matlabbatch{1}.spm.spatial.preproc.warp.mrf     = 1;             
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg     = 4;
matlabbatch{1}.spm.spatial.preproc.warp.affreg  = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm    = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp    = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write   = [0 0];
spm_jobman('run', matlabbatch);

% Build single-scan mask and stripped T1
c1_1 = fullfile(outdir, ['c1' spm_file(t1_path_1,'filename')]);       % GM probability map
c2_1 = fullfile(outdir, ['c2' spm_file(t1_path_1,'filename')]);       % WM prob map
c3_1 = fullfile(outdir, ['c3' spm_file(t1_path_1,'filename')]);       % CSF prob map
gm1  = spm_read_vols(spm_vol(c1_1));
wm1  = spm_read_vols(spm_vol(c2_1));
csf1 = spm_read_vols(spm_vol(c3_1));
mask_single = (gm1 + wm1 + csf1) > 0.3;                               % Adding the prob maps gives an estimate of the voxel’s brain probability
                                                                      % the 0.3 threshold keeps the model where it's believed that there's at 
                                                                      % least 30% prob of tissues being brain tissues 

orig_vol = spm_vol(m1_path);                                          % Use bias-corrected T1 (m*) before masking.
orig_dat = spm_read_vols(orig_vol);
single_stripped_dat      = orig_dat .* mask_single;                   % element-wise mask
orig_vol.fname           = single_stripped_fn;
spm_write_vol(orig_vol, single_stripped_dat);

% B) Mean of the two T1s, then segment the mean -> c1/c2/c3 of mean
V1 = spm_vol(t1_path_1); I1 = spm_read_vols(V1);
V2 = spm_vol(t1_path_2); I2 = spm_read_vols(V2);
Imean = (I1 + I2) / 2;
Vmean = V1; Vmean.fname = mean_t1_path; spm_write_vol(Vmean, Imean);

matlabbatch = [];
matlabbatch{1}.spm.spatial.preproc.channel.vols     = {[Vmean.fname ',1']};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg  = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write    = [0 1];
for t = 1:3
    matlabbatch{1}.spm.spatial.preproc.tissue(t).tpm    = {fullfile(spm('Dir'),'tpm',['TPM.nii,' num2str(t)])};
    matlabbatch{1}.spm.spatial.preproc.tissue(t).ngaus  = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(t).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(t).warped = [0 0];
end
for t = 4:6
    matlabbatch{1}.spm.spatial.preproc.tissue(t).tpm    = {fullfile(spm('Dir'),'tpm',['TPM.nii,' num2str(t)])};
    matlabbatch{1}.spm.spatial.preproc.tissue(t).ngaus  = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(t).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(t).warped = [0 0];
end
matlabbatch{1}.spm.spatial.preproc.warp.mrf     = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg     = 4;
matlabbatch{1}.spm.spatial.preproc.warp.affreg  = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm    = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp    = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write   = [0 0];
spm_jobman('run', matlabbatch);

% Build mask from mean T1 tissue maps
gm = spm_read_vols(spm_vol(fullfile(outdir,'c1mean_T1.nii')));
wm = spm_read_vols(spm_vol(fullfile(outdir,'c2mean_T1.nii')));
csf= spm_read_vols(spm_vol(fullfile(outdir,'c3mean_T1.nii')));
mask = (gm + wm + csf) > 0.3;

mask_vol       = spm_vol(fullfile(outdir,'c1mean_T1.nii'));
mask_vol.fname = brain_mask_path;
spm_write_vol(mask_vol, mask);

% Apply mask to bias-corrected mean to produce the final brain-only mean T1
stripped_vol       = spm_vol(fullfile(outdir,'mmean_T1.nii'));
stripped_data      = spm_read_vols(stripped_vol) .* mask;
stripped_vol.fname = mean_t1_stripped;
spm_write_vol(stripped_vol, stripped_data);

% -------- QC figures (normal mode)

% 2x2 comparison
try
    mean_orig        = spm_read_vols(spm_vol(mean_t1_path));
    orig_t1_biascorr = spm_read_vols(spm_vol(m1_path));
    single_stripped  = single_stripped_dat;
    sl               = round(size(orig_t1_biascorr,3)/2);

    f2 = figure('Name','Comparison of Original vs Skull-Stripped T1','Color','w');
    sgtitle('Comparison of Original vs Skull-Stripped T1');
    subplot(2,2,1); imshow(mean_orig(:,:,sl), []);        title('Mean T1');
    subplot(2,2,2); imshow(stripped_data(:,:,sl), []);    title('Mean T1 (stripped)');
    subplot(2,2,3); imshow(orig_t1_biascorr(:,:,sl), []); title('T1#1 (bias-corrected)');
    subplot(2,2,4); imshow(single_stripped(:,:,sl), []);  title('T1#1 (stripped)');
    
    if opt.SaveFigures && ~isempty(opt.SaveDir)
        png_out = fullfile(opt.SaveDir,[prefix 'T1_comparison.png']);
        fig_out = fullfile(opt.SaveDir,[prefix 'T1_comparison.fig']);
        exportgraphics(f2, png_out, 'Resolution', 200);
        savefig(f2, fig_out);
        fprintf('Saved in:\n  %s\n  %s\n', png_out, fig_out);
    end
catch ME
    warning(ME.identifier,'Step 2 comparison figure not saved: %s', ME.message);
end
end

