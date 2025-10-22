function resliced_mask_path = reslice_t1mask_to_epi(mean_epi_path, brain_mask_path) 
% Inputs:
%   mean_epi_path    (char) — path to mean EPI (reference image)
%   brain_mask_path  (char) — path to T1-space brain mask (source image, label map)
%
% Output:
%   resliced_mask_path (char) — path to mask resliced onto the mean-EPI grid, prefixed with 'r'
%
% Purpose:
%   Resample a T1-space intracranial mask onto the mean-EPI grid using SPM’s
%   coreg.write. This uses nearest-neighbor (NN) interpolation to preserve
%   discrete labels.
%
%
% Author: Nagham Nessim
% University of Geneva, 2025

% Reslice a T1-space brain mask onto the mean-EPI grid (NN), writes 'r*' file.

    mean_epi_path   = char(mean_epi_path);
    brain_mask_path = char(brain_mask_path);

    [pm,fm,ex] = fileparts(brain_mask_path);
    resliced_mask_path = fullfile(pm, ['r' fm ex]);

    if ~isfile(resliced_mask_path)
        spm('defaults','FMRI'); spm_jobman('initcfg');
        matlabbatch = [];
        matlabbatch{1}.spm.spatial.coreg.write.ref                 = { [mean_epi_path ',1'] };     % Reference grid = mean EPI. Write onto this grid
        matlabbatch{1}.spm.spatial.coreg.write.source              = { [brain_mask_path ',1'] };   % Source = T1-space mask (to be resampled)

        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp     = 0;                            % The brain mask is a label/categorical image (0/1 or multi-label). 
                                                                                                   % Linear, spline, or other interpolation would create non-integer 
                                                                                                   % hybrids and shrink/blur boundaries. Nearest-neighbor preserves
                                                                                                   % labels exactly and keeps the mask binary. Nearest-neighbor (NN) 
                                                                                                   % picks the single closest voxel’s value.
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap       = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.write.roptions.mask       = 0;
        matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix     = 'r';                          % Prefix for output files
        spm_jobman('run', matlabbatch);
    end

end
