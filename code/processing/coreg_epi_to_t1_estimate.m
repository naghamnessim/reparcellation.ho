function realigned_path_out = coreg_epi_to_t1_estimate(mean_t1_stripped, mean_epi_path, realigned_path)
% Inputs:
%   mean_t1_stripped (char) — skull-stripped T1 (reference image, fixed space)
%   mean_epi_path    (char) — mean EPI (source image, to be aligned)
%   realigned_path   (char) — 4D motion-corrected EPI (other images whose headers will be updated)
%
% Output:
%   realigned_path_out (char) — same as input realigned_path; data unchanged, headers updated in-place
%
% Purpose:
%   Perform header-only rigid coregistration of the EPI to the skull-stripped T1
%   using normalized mutual information (NMI). The estimated rigid transform is applied to the
%   headers of both the mean EPI and the whole 4D EPI series, avoiding an extra interpolation now
%   so reslicing happens only once later in the pipeline.
%
% Author: Nagham Nessim
% University of Geneva, 2025


    % Ensure char (SPM prefers char arrays)
    mean_t1_stripped = char(mean_t1_stripped);
    mean_epi_path    = char(mean_epi_path);
    realigned_path   = char(realigned_path);

    % Minimal checksto fail fast if inputs are missing
    assert(isfile(mean_t1_stripped), 'Missing T1: %s', mean_t1_stripped);
    assert(isfile(mean_epi_path),    'Missing mean EPI: %s', mean_epi_path);
    assert(isfile(realigned_path),   'Missing 4D EPI: %s', realigned_path);

    % Build and run SPM batch (estimate-only; updates headers) 
    % “estimate” computes 6-DOF rigid transform that maps the source (mean EPI) into the reference (skull-stripped T1) using an NMI cost function.
    % Header-only coreg keeps voxel data for a single reslice (one interpolation).
    % no new intensities are computed, no interpolation is applied. The files now point to a new pose in world space, but the raw voxel grid and values are untouched.

    spm('defaults','FMRI'); spm_jobman('initcfg');
    matlabbatch = [];
    matlabbatch{1}.spm.spatial.coreg.estimate.ref    = { [mean_t1_stripped ',1'] };       % Reference (fixed) = T1 brain
    matlabbatch{1}.spm.spatial.coreg.estimate.source = { [mean_epi_path ',1'] };          % Source (moving) = mean EPI
    matlabbatch{1}.spm.spatial.coreg.estimate.other  = { realigned_path };                % Apply header transform to all frames of 4D EPI
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol      = ...
        [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];                  % Gaussian smoothing
    spm_jobman('run', matlabbatch);

    % Return same path (data unchanged, headers updated)
    realigned_path_out = realigned_path;
end

