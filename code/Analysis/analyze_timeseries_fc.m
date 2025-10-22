function [roi_timeseries, roi_ts_z, fc_matrix, roi_labels] = analyze_timeseries_fc( ...
    warped_epi_path, ho_resliced_path, varargin)

% Inputs:
%   warped_epi_path  (char) — path to 4D EPI in MNI space, shape [X Y Z T]
%   ho_resliced_path (char) — path to HO atlas resliced to the EPI grid
%   Name-Value:
%     'SaveFigures'  (logical, default false) — save PNGs
%     'SaveDir'      (char, default '')       — output folder for figures
%     'FigurePrefix' (char, default '')       — prefix for figure filenames
%
% Outputs:
%   roi_timeseries ([T x R] double) — raw mean time series per ROI (columns are ROIs)
%   roi_ts_z      ([T x R] double) — per-ROI z-scored time series across time (mean 0, std 1)
%   fc_matrix     ([R x R] double) — Pearson correlation between ROI time series
%   roi_labels    ([R x 1] double) — atlas labels corresponding to ROI
%
% Purpose:
%   Extract parcel-wise time series from an atlas-aligned fMRI (MNI space) by averaging
%   voxel signals within each ROI, then compute an ROI×ROI Pearson correlation matrix
%   as a simple functional connectivity (FC) estimate.
% Author: Nagham Nessim
% University of Geneva, 2025

p = inputParser;
addParameter(p,'SaveFigures',false,@islogical);
addParameter(p,'SaveDir','',@(s)ischar(s)||isstring(s));
addParameter(p,'FigurePrefix','',@(s)ischar(s)||isstring(s));
parse(p,varargin{:});
opt = p.Results;

if ~isempty(opt.SaveDir) && ~exist(opt.SaveDir,'dir')
    mkdir(opt.SaveDir);
end
prefix = char(opt.FigurePrefix);
if ~isempty(prefix) && ~endsWith(prefix,'_')
    prefix = [prefix '_'];
end

% Load Data
V_func = spm_vol(warped_epi_path);
func_data = spm_read_vols(V_func);                                          % Read all voulmes of the 4D EPI [x, y, z, t]

V_atlas = spm_vol(ho_resliced_path);
atlas_data = spm_read_vols(V_atlas);                                        % Label 3D image with integer ROI IDs [x, y, z]


% Safety Check: confirm grid match
assert(isequal(size(func_data,1),size(atlas_data,1)) && ...
       isequal(size(func_data,2),size(atlas_data,2)) && ...
       isequal(size(func_data,3),size(atlas_data,3)), ...
       'Atlas/EPI grid mismatch — ensure HO was resliced onto EPI grid correctly.'); 
disp("Func data size:");  disp(size(func_data));
disp("Atlas data size:"); disp(size(atlas_data));

% ROI labels
roi_labels = unique(atlas_data);
roi_labels(roi_labels == 0 | isnan(roi_labels)) = [];                       % remove background/NaN (HO uses 0 for outside-atlas)
n_rois = numel(roi_labels);
n_timepoints = size(func_data, 4);

% extract ROI-wise mean time series
roi_timeseries = zeros(n_timepoints, n_rois);
for r = 1:n_rois
    roi_mask = (atlas_data == roi_labels(r));
    for t = 1:n_timepoints
        vol_t = func_data(:, :, :, t);
        roi_timeseries(t, r) = mean(vol_t(roi_mask), 'omitnan');
    end
end

% z-score per ROI across time
roi_ts_z = zscore(roi_timeseries, 0, 1);                                    % Standardize each ROI’s time series (mean 0, std 1) 
                                                                            % comparable scales prior to Functional Connectivity; Pearson coefficient r is 
                                                                            % scale-invariant but z-scoring aids visualization and numerical stability

% ---------- Figure 1: overview ----------
f1 = figure('Name', 'ROI Time Series (MNI Space)', 'Color', 'w');
subplot(2,1,1)
plot(roi_timeseries); 
title('Raw ROI Time Series'); 
xlabel('Time (TR)'); 
ylabel('Signal')

subplot(2,1,2)
plot(roi_ts_z); 
title('Z-scored ROI Time Series'); 
xlabel('Time (TR)'); 
ylabel('Z-score')

if opt.SaveFigures && ~isempty(opt.SaveDir)
    png_out = fullfile(opt.SaveDir, sprintf('%sroi_timeseries_overview.png', prefix));
    fig_out = fullfile(opt.SaveDir, sprintf('%sroi_timeseries_overview.fig', prefix));
    exportgraphics(f1, png_out, 'Resolution', 200);
    savefig(f1, fig_out);
    fprintf('Raw ROI Timeseries Saved in : %s\n          %s\n', png_out, fig_out);
end

% ---------- Figure 2: z-scored per ROI ----------
n_cols = 8;
n_rows = ceil(n_rois / n_cols);
f2 = figure('Name', 'Z-scored ROI Time Series', 'Color', 'w', 'Position', [100, 100, 1600, 900]);
for r = 1:n_rois
    subplot(n_rows, n_cols, r)
    plot(roi_ts_z(:, r), 'b', 'LineWidth', 1.2)
    xlim([1 n_timepoints])
    xt = unique(round(linspace(1, n_timepoints, min(5, n_timepoints))));
    set(gca, 'XTick', xt)
    title(sprintf('ROI #%d', r))
    grid on
end
sgtitle('Z-scored Time Series per ROI')
if opt.SaveFigures && ~isempty(opt.SaveDir)
    png_out = fullfile(opt.SaveDir, sprintf('%sroi_timeseries_z_byROI.png', prefix));
    fig_out = fullfile(opt.SaveDir, sprintf('%sroi_timeseries_z_byROI.fig', prefix));
    exportgraphics(f2, png_out, 'Resolution', 200);
    savefig(f2, fig_out);
    fprintf('Z-Scored Timeseries Saved in: %s\n          %s\n', png_out, fig_out);
end

% ---------- Figure 3: raw per ROI ----------
f3 = figure('Name', 'Raw ROI Time Series', 'Color', 'w', 'Position', [100, 100, 1600, 900]);
for r = 1:n_rois
    subplot(n_rows, n_cols, r)
    plot(roi_timeseries(:, r), 'b', 'LineWidth', 1.2)
    xlim([1 n_timepoints])
    xt = unique(round(linspace(1, n_timepoints, min(5, n_timepoints))));
    set(gca, 'XTick', xt)
    title(sprintf('ROI #%d', r))
    grid on
end
sgtitle('Raw Time Series per ROI')
if opt.SaveFigures && ~isempty(opt.SaveDir)
    png_out = fullfile(opt.SaveDir, sprintf('%sroi_timeseries_raw_byROI.png', prefix)); 
    fig_out = fullfile(opt.SaveDir, sprintf('%sroi_timeseries_raw_byROI.fig', prefix));
    exportgraphics(f3, png_out, 'Resolution', 200);
    savefig(f3, fig_out);
    fprintf('Raw Timeseries per ROI Saved in: %s\n          %s\n', png_out, fig_out);
end

% ---------- FC Calculation----------
fc_matrix = corr(roi_ts_z);                                                 % Pearson correlation among ROI columns. matrix size wil be [n_rois x n_rois]

% ---------- Figure 4: FC matrix ----------
f4 = figure('Name', 'Functional Connectivity Matrix', 'Color', 'w');
imagesc(fc_matrix); axis square
colorbar; colormap(jet)
title('Functional Connectivity')
xlabel('ROI'); ylabel('ROI');
if opt.SaveFigures && ~isempty(opt.SaveDir)
    png_out = fullfile(opt.SaveDir, sprintf('%sfc_matrix.png', prefix));
    fig_out = fullfile(opt.SaveDir, sprintf('%sfc_matrix.fig', prefix));
    exportgraphics(f4, png_out, 'Resolution', 200);
    savefig(f4, fig_out);
    fprintf('FC Matrix Saved in: %s\n          %s\n', png_out, fig_out);
end
end

