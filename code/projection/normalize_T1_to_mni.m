function [y_field_path, iy_field_path] = normalize_T1_to_mni(mean_t1_stripped_path)
% Inputs:
%   mean_t1_stripped_path (char) — bias-corrected, skull-stripped mean T1 (native space)
%
% Outputs:
%   y_field_path  (char) — forward deformation field  (native → MNI), file: y_<name>.nii
%   iy_field_path (char) — inverse deformation field  (MNI → native),  file: iy_<name>.nii
%
% Purpose:
%   Run SPM’s (unified segmentation + normalization) on the subject’s T1 to estimate a
%   nonlinear warp between native space and MNI152 template space. Writes both the
%   forward (y_) and inverse (iy_) fields so it's possible to move images to MNI space or bring atlas
%   labels into subject space later with a single resampling.
%

spm('defaults', 'FMRI');
spm_jobman('initcfg');

% Prepare batch
matlabbatch = [];
matlabbatch{1}.spm.spatial.preproc.channel.vols = {[mean_t1_stripped_path ',1']};      % Input channel = T1; ',1' selects the 1st and only volume.
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;                            % bias-field regularization
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;                              % bias smoothness
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];                              % save the bias-corrected image (*m)

for t = 1:6                                                                            % Tissue posterior prob maps (GM, WM, CSF, bone, soft tissue, air)
    matlabbatch{1}.spm.spatial.preproc.tissue(t).tpm = {fullfile(spm('Dir'), 'tpm', ['TPM.nii,' num2str(t)])};
    matlabbatch{1}.spm.spatial.preproc.tissue(t).ngaus = 1;                            % 1 gaussian per tissue class
    matlabbatch{1}.spm.spatial.preproc.tissue(t).native = [1 0];                       % write native space for tissue maps to build mask
    matlabbatch{1}.spm.spatial.preproc.tissue(t).warped = [0 0];
end

matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = 4;
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];                                 % Write both forward y_* (native to MNI) 
                                                                                       % and inverse iy_* (bring atlas labels to native subject space) fields.

% Run segmentation
spm_jobman('run', matlabbatch);

% Output paths
[p, name, ~] = fileparts(mean_t1_stripped_path);
y_field_path = fullfile(p, ['y_' name '.nii']);                                        % Forward deformation field file
iy_field_path = fullfile(p, ['iy_' name '.nii']);                                      % Inverse deformation field file

end
