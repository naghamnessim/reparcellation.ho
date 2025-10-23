function [warped_epi_path, rT1_on_wEPI] = warp_epi_to_mni(epi_stripped_path, mean_t1_stripped, y_field_path) 
% Inputs:
%   epi_stripped_path (char) — skull-stripped 4D EPI in native space (analysis EPI)
%   mean_t1_stripped  (char) — bias-corrected, skull-stripped mean T1 (native space)
%   y_field_path      (char) — forward deformation field y_* (native → MNI) from unified segmentation
%
% Outputs:
%   warped_epi_path (char) — EPI warped to MNI space (prefixed 'w', 2 mm voxels)
%   rT1_on_wEPI     (char) — T1 warped to MNI then resliced onto the wEPI grid (prefixed 'rEPI_')
%
% Purpose:
%   (1) Apply the **forward nonlinear warp** (native→MNI) to the full 4D EPI and write at 2 mm,
%   (2) Warp the T1 to MNI and **reslice it onto the exact wEPI sampling grid** for 1:1 overlays.
%
%
% Author: Nagham Nessim
% University of Genova, 2025


    % Normalize inputs & check for files' existence
    epi_stripped_path = char(epi_stripped_path);
    mean_t1_stripped  = char(mean_t1_stripped);
    y_field_path      = char(y_field_path);

    assert(isfile(epi_stripped_path), 'EPI not found: %s', epi_stripped_path);
    assert(isfile(mean_t1_stripped),  'T1 not found: %s', mean_t1_stripped);
    assert(isfile(y_field_path),      'Deformation field not found: %s', y_field_path);

    % SPM setup 
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');

    %-------------------
    % Warp the entire 4D EPI into MNI space (2 mm)
    % ------------------
    % Build a list of all timepoints so the whole time series is warped.
    Vepi = spm_vol(epi_stripped_path);        % one header per timepoint
    nF   = numel(Vepi);
    resample_epi = cell(nF,1);
    for i = 1:nF
        resample_epi{i} = sprintf('%s,%d', epi_stripped_path, i);
    end

    matlabbatch = [];
    matlabbatch{1}.spm.spatial.normalise.write.subj.def        = {y_field_path};            % Use forward deformation y_* (subject native to MNI space)
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample   = resample_epi;              % Resample the full 4D EPI
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb     = [-78 -112 -70; 78 76 85];  % MNI152 bounding box
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox    = [2 2 2];                   % match common atlas resolution (2 mm isotropic)
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;                         % 4th-degree B-spline interpolation for continuous EPI intensities

    spm_jobman('run', matlabbatch);

    [pE, nameE, ~]  = fileparts(epi_stripped_path); 
    warped_epi_path = fullfile(pE, ['w' nameE '.nii']);                                     % prefix 'w' for output file (normalized to template) 
    assert(isfile(warped_epi_path), 'Warped EPI not created: %s', warped_epi_path);

    % -----------------
    % Warp T1 to MNI (2 mm), then reslice T1 onto the warped EPI grid
    % -----------------
    % a) T1 to MNI using the same forward field (Subject to MNI)
    matlabbatch = [];
    matlabbatch{1}.spm.spatial.normalise.write.subj.def        = {y_field_path};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample   = {mean_t1_stripped};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb     = [-78 -112 -70; 78 76 85];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox    = [2 2 2];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;

    spm_jobman('run', matlabbatch);

    [pT, nameT, ~] = fileparts(mean_t1_stripped);
    warped_t1_path = fullfile(pT, ['w' nameT '.nii']);
    assert(isfile(warped_t1_path), 'Warped T1 not created: %s', warped_t1_path);

    % b) reslice warped T1 onto the exact sampling grid of warped EPI (wEPI) (no new warp)
    %Use vol-1 of wEPI as the reference grid.
    matlabbatch = [];
    matlabbatch{1}.spm.spatial.coreg.write.ref                 = {[warped_epi_path ',1']};    % Reference grid = wEPI (target sampling
    matlabbatch{1}.spm.spatial.coreg.write.source              = {[warped_t1_path ',1']};     % Source = warped T1 (already in MNI space); just reslice to wEPI grid.
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp     = 4;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap       = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask       = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix     = 'rEPI_';

    spm_jobman('run', matlabbatch);

    rT1_on_wEPI = fullfile(pT, ['rEPI_w' nameT '.nii']);  % Matches wEPI affine matrix exactly. Perfect for overlays & voxelwise compare.
    assert(isfile(rT1_on_wEPI), 'Resliced T1 (on wEPI grid) not created: %s', rT1_on_wEPI);
end

