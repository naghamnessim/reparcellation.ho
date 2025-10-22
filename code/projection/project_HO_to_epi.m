function ho_resliced_path = project_HO_to_epi(warped_epi_path, ho_gz)
% Inputs:
%   warped_epi_path (char) — path to MNI-warped functional image (w*.nii). If 4D, vol 1 is used as reference.
%   ho_gz          (char) — path to HarvardOxford-cort-maxprob-thr25-2mm.nii(.gz)
%
% Output:
%   ho_resliced_path (char) — path to HO atlas resliced onto the exact wEPI grid (prefixed 'r')
%
% Purpose:
%   Reslice (grid match only) the Harvard–Oxford atlas onto the MNI-warped EPI grid using
%   SPM’s Coreg:Write. Uses nearest-neighbor interpolation to preserve parcel labels.
%
% Author: Nagham Nessim
% University of Geneva, 2025
   
    ho_atlas_path  = ho_gz;
    assert(isfile(ho_atlas_path), 'HO atlas not found: %s', ho_atlas_path);

    % SPM setup
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');

    % Use vol 1 of wEPI as reference grid
    ref = sprintf('%s,1', warped_epi_path);                                                  % Reference = wEPI grid
    src = sprintf('%s,1', ho_atlas_path);                                                    % Source = HO atlas (label map) to be put on the wEPI sampling grid

    % Reslice HO to wEPI grid (Nearest Neighbors to preserve labels (binary)
    matlabbatch = [];
    matlabbatch{1}.spm.spatial.coreg.write.ref                 = {ref};
    matlabbatch{1}.spm.spatial.coreg.write.source              = {src};
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp     = 0;                          % NN interpolation
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap       = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask       = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix     = 'r';                        % prefix for the output file

    spm_jobman('run', matlabbatch);

    % Output name, and check sizes match wEPI
    [p, name, ~]    = fileparts(ho_atlas_path);
    ho_resliced_path = fullfile(p, ['r' name '.nii']);
    assert(isfile(ho_resliced_path), 'Resliced HO not created: %s', ho_resliced_path);

    Vref = spm_vol(warped_epi_path); if numel(Vref)>1, Vref = Vref(1); end
    Vho  = spm_vol(ho_resliced_path);
    szE  = Vref.dim;                                                                        % [X Y Z] of wEPI grid.
    szH  = Vho.dim;                                                                         % [X Y Z] of resliced HO grid.
    assert(isequal(szE, szH), 'Grid mismatch after reslice: EPI [%s] vs HO [%s].', ...      % check sizes match
           num2str(szE), num2str(szH));
end

