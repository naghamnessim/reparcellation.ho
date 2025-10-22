function [realigned_path, mean_epi_path, rp_path] = preprocess_functional(func_nii, json_path, varargin)
% Inputs:
%   func_nii  (char) — full path to the functional NIfTI (3D/4D), native EPI space
%   json_path (char) — full path to the BIDS JSON sidecar for this run (holds TR, SliceTiming)
%
% Outputs:
%   realigned_path (char) — path to motion-corrected 4D NIfTI (prefix 'ra')
%   mean_epi_path  (char) — path to the mean EPI produced during realignment (prefix 'meana')
%
% Purpose:
%   Perform slice-timing correction (STC) followed by rigid-body motion correction (realignment)
%   in SPM for a 4D fMRI time series, and return the key derivative paths for downstream steps.
%

p = inputParser;
addParameter(p,'SaveFigures',false,@islogical);
addParameter(p,'SaveDir','',@(s)ischar(s)||isstring(s));
addParameter(p,'FigurePrefix','step1_',@(s)ischar(s)||isstring(s));
addParameter(p,'QCOnly',false,@islogical);
addParameter(p,'CachedRealignedPath','',@(s)ischar(s)||isstring(s));
parse(p,varargin{:});
opt = p.Results;

if opt.SaveFigures && ~isempty(opt.SaveDir) && ~exist(opt.SaveDir,'dir')
    mkdir(opt.SaveDir);
end
prefix = char(opt.FigurePrefix);
if ~isempty(prefix) && ~endsWith(prefix,'_'), prefix = [prefix '_']; end

% -----------------------
% QC-only mode: only regenerate motion plots from cached outputs
% -----------------------
if opt.QCOnly
    if isempty(opt.CachedRealignedPath) || ~isfile(opt.CachedRealignedPath)
        error('QCOnly requested but path is missing or invalid.');
    end

    % reconstruct paths from cached outputs
    realigned_path = opt.CachedRealignedPath;
    [out_dir, ra_name, ~] = fileparts(realigned_path);
    if startsWith(ra_name,'ra')
        base_a = ra_name(2:end);
    elseif startsWith(ra_name,'r')
        base_a = ['a' ra_name(2:end)];
    else
        base_a = ['a' ra_name];
    end

    rp_path = fullfile(out_dir, ['rp_' base_a '.txt']);
    mean_epi_path = fullfile(out_dir, ['meana' spm_file(func_nii, 'filename')]);

    % read TR from JSON
    fid = fopen(json_path); raw = fread(fid, inf); str = char(raw'); fclose(fid);
    meta = jsondecode(str);
    TR = meta.RepetitionTime;

    % Motion Parameter visualization only
    if opt.SaveFigures && isfile(rp_path)
        try
            mp = readmatrix(rp_path);
            t = (0:size(mp,1)-1).' * TR;
            trans   = mp(:,1:3);
            rot_deg = mp(:,4:6) * (180/pi);

            fig = figure('Name','SPM Realignment Parameters','Color','w',...
                         'Units','pixels','Position',[100 100 1100 600]);
            tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

            nexttile;
            plot(t, trans, 'LineWidth', 1.5); grid on;
            xlabel('Time (s)'); ylabel('Translation (mm)');
            legend({'X','Y','Z'}, 'Location','best'); title('Head Translations');

            nexttile;
            plot(t, rot_deg, 'LineWidth', 1.5); grid on;
            xlabel('Time (s)'); ylabel('Rotation (deg)');
            legend({'Pitch (x)','Roll (y)','Yaw (z)'}, 'Location','best'); title('Head Rotations');

            png_out = fullfile(opt.SaveDir, [prefix 'motion_params''.png']);
            exportgraphics(fig, png_out, 'Resolution', 150);
            savefig(fig, fullfile(opt.SaveDir, [prefix 'motion_params''.fig']));
            fprintf('[QC-only] Saved motion plots: %s\n', png_out);
        catch ME
            warning(ME.identifier,'QC-only visualization failed: %s', ME.message);
        end
    end

    return;
end

% -----------------------
% BIDS JSON LOADING
% -----------------------
if ~isfile(json_path), error('BIDS JSON not found: %s', json_path); end 
fid = fopen(json_path);
raw = fread(fid, inf);
str = char(raw');
fclose(fid);
meta = jsondecode(str);

% -----------------------
% DERIVE Slice Time Correction PARAMETERS
% -----------------------
TR = meta.RepetitionTime;
slice_timing = meta.SliceTiming;
nslices = length(slice_timing);
TA = TR - TR / nslices;                                                              % SPM temporal acquisition window (TR minus one slice interval)
[~, slice_order] = sort(slice_timing);
refslice = round(nslices / 2);                                                       % Extract Middle slice

% -----------------------
% Slice timing correction btach
% -----------------------
hdr = spm_vol(func_nii);                                                             % Load volume headers to enumerate 4D timepoints.
n_vols = length(hdr);                                                                % Number of timepoints (T) in series.
scans = cell(n_vols, 1);                                                            
for i = 1:n_vols
    scans{i} = [func_nii ',' num2str(i)];                                            
end

matlabbatch{1}.spm.temporal.st.scans = {scans};                                      % Provide all volumes for STC.
matlabbatch{1}.spm.temporal.st.nslices = nslices;                                    % Defines slice grid for resampling.
matlabbatch{1}.spm.temporal.st.tr = TR;                                              % TR informs timing grid between volumes.
matlabbatch{1}.spm.temporal.st.ta = TA;                                              % Acquisition window
matlabbatch{1}.spm.temporal.st.so = slice_order;                                     % Slice acquisition order (ascending)
matlabbatch{1}.spm.temporal.st.refslice = refslice;                                  % Align to middle slice time
matlabbatch{1}.spm.temporal.st.prefix = 'a';                                         % 'a' prefix marks slice-timed images.

spm('defaults', 'FMRI');                                                              
spm_jobman('initcfg');                                                               
spm_jobman('run', matlabbatch);                                                      

% -----------------------
% Realignment (Motion Correction)
% Estimate 6-DOF motion and reslice series
% -----------------------
corrected_path = fullfile(fileparts(func_nii), ['a' spm_file(func_nii, 'filename')]);  %time-corrected path
hdr_corr = spm_vol(corrected_path);
scans = cell(numel(hdr_corr), 1);
for i = 1:numel(hdr_corr)
    scans{i} = [corrected_path ',' num2str(i)];
end

clear matlabbatch;                                                                    % Avoid cross-talk with previous batch fields.
matlabbatch{1}.spm.spatial.realign.estwrite.data = {scans};                           % “estwrite” estimates motion and reslices.
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.99;                  % High sampling of cost function — more accurate, slower.
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 3;                         % Sampling separation (mm)
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 6;                        % Smoothing (mm) on images during estimation.
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;                         % Register to first image (0) vs mean (1)
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;                      % 2nd-degree B-spline interpolation during parameter estimation.
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];                   % Write resliced all images (2) and mean image (1).
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;                      % 4th-degree B-spline interpolation for final reslicing.
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];                  % No wrap
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;                        % Mask edges
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';                    % 'r' prefix marks motion-corrected outputs.

spm_jobman('run', matlabbatch);                                                       

% Output paths
realigned_path = fullfile(fileparts(func_nii), ['ra' spm_file(func_nii, 'filename')]);
mean_epi_path = fullfile(fileparts(func_nii), ['meana' spm_file(func_nii, 'filename')]);

% SPM writes motion parameters for the stack it realigned
stc_basename = spm_file(corrected_path, 'basename');
rp_path      = fullfile(fileparts(func_nii), ['rp_' stc_basename '.txt']); 

% Visulaization of Motion Parameters (translation and rotation)
if opt.SaveFigures && isfile(rp_path)
    try
        mp = readmatrix(rp_path);
        t = (0:size(mp,1)-1).' * TR;
        trans   = mp(:,1:3);
        rot_deg = mp(:,4:6) * (180/pi);

        fig = figure('Name','SPM Realignment Parameters','Color','w',...
                     'Units','pixels','Position',[100 100 1100 600]);
        tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

        nexttile;
        plot(t, trans, 'LineWidth', 1.5); grid on;
        xlabel('Time (s)'); ylabel('Translation (mm)');
        legend({'X','Y','Z'}, 'Location','best'); title('Head Translations');

        nexttile;
        plot(t, rot_deg, 'LineWidth', 1.5); grid on;
        xlabel('Time (s)'); ylabel('Rotation (deg)');
        legend({'Pitch (x)','Roll (y)','Yaw (z)'}, 'Location','best'); title('Head Rotations');

        png_out = fullfile(opt.SaveDir, [prefix 'motion_params''.png']);
        exportgraphics(fig, png_out, 'Resolution', 150);
        savefig(fig, fullfile(opt.SaveDir, [prefix 'motion_params''.fig']));
        fprintf('Saved: %s\n', png_out);
    catch ME
        warning(ME.identifier,'Motion visualization failed: %s', ME.message);
    end
end

end
