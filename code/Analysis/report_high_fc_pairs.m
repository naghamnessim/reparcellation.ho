function [roi_pair_table] = report_high_fc_pairs( ...
    fc_matrix, ho_resliced_path, base_path, k_std, varargin)
% Inputs:
%   fc_matrix        ([R x R] double) — symmetric Pearson-r FC matrix between ROIs
%   ho_resliced_path (char)           — path to HO atlas on EPI grid (used to recover ROI IDs present)
%   base_path        (char)           — project root; used to look up HO XML label file
%   k_std            (double)         — threshold in SD (standard deviation) units above mean (mean + k*std)
%   Name-Value:
%     'SaveFigures'  (logical, default false) — (NOTE: current implementation does not draw a figure)
%     'SaveTable'    (logical, default false) — write CSV of high-FC pairs
%     'SaveDir'      (char/string, default '') — output folder for figure/CSV
%     'FigurePrefix' (char/string, default '') — filename prefix
%
% Outputs:
%   roi_pair_table (table) — columns: ROI_1, Label_1, ROI_2, Label_2, Correlation
%
% Purpose:
%   Identify highly correlated ROI pairs using a simple threshold:
%   r ≥ mean(r_offdiag) + k_std·std(r_offdiag). Maps FC indices to HO atlas labels
%   (via XML) and returns a table, optionally saved as CSV.

p = inputParser;
addParameter(p,'SaveFigures',false,@islogical);
addParameter(p,'SaveTable',false,@islogical);
addParameter(p,'SaveDir','',@(s)ischar(s)||isstring(s));
addParameter(p,'FigurePrefix','',@(s)ischar(s)||isstring(s));
parse(p,varargin{:});
opt = p.Results;

if ~isempty(opt.SaveDir) && ~exist(opt.SaveDir,'dir')
    mkdir(opt.SaveDir);
end
prefix = char(opt.FigurePrefix);
if ~isempty(prefix) && ~endsWith(prefix,'_'), prefix = [prefix '_']; 
end

% read HO labels (XML)
% Default location under project root
xml_path = fullfile(base_path, 'Harvard Oxford Atlas', 'HarvardOxford-Cortical.xml');
if ~isfile(xml_path)
    xml_path = fullfile(fileparts(ho_resliced_path), 'HarvardOxford-Cortical.xml');   % if it does not exist, fallback to the same folder as the resliced HO atlas

if isfile(xml_path)
    try
        xDoc = xmlread(xml_path);
        labels = xDoc.getElementsByTagName('label');
        ho_labels = cell(labels.getLength, 1);
        for i = 0:labels.getLength-1
            ho_labels{i+1} = char(labels.item(i).getFirstChild.getData);
        end
    catch ME
        warning(ME.identifier, 'Could not parse HO XML; using generic labels. %s', ME.message);
        ho_labels = {};
    end
else
    warning('report_high_fc_pairs:NoXML', 'HO XML not found (%s). Using generic labels.', xml_path);
    ho_labels = {};
end

% -------------------- atlas ROI ids (to map FC indices to HO IDs) --------------------
atlas_img = spm_read_vols(spm_vol(ho_resliced_path));
roi_ids = unique(atlas_img(:));
roi_ids(roi_ids == 0 | isnan(roi_ids)) = [];                                                 % Remove background and NaNs.
roi_ids = sort(roi_ids(:));                                                                  % column, ascending order

% Check that FC size matches the number of ROIs derived from atlas image
if size(fc_matrix,1) ~= numel(roi_ids) || size(fc_matrix,2) ~= numel(roi_ids)
    warning('report_high_fc_pairs:SizeMismatch', ...
        'FC is %dx%d but atlas has %d ROIs; label mapping may be wrong.', ...
        size(fc_matrix,1), size(fc_matrix,2), numel(roi_ids));                               % Mismatch suggests different ROI sets
end

% -------------------- Apply threshold & pairs --------------------
fc_no_diag = fc_matrix;                                                                      % copy the matrix to preserve the original one
fc_no_diag(eye(size(fc_matrix)) == 1) = NaN;                                                 % exclude self-correlations

fc_vals = fc_no_diag(~isnan(fc_no_diag));
threshold = mean(fc_vals, 'omitnan') + k_std * std(fc_vals, 'omitnan');                      % threshold r ≥ mean + k*std

[row_idx, col_idx] = find(fc_no_diag >= threshold);                                          % Get all supra-threshold pairs.
n_pairs = numel(row_idx);                                                                    % number of entries meeting the threshold

% -------------------- build table --------------------
roi_pair_matrix = cell(n_pairs, 5);                                                          % ROI_1, Label_1, ROI_2, Label_2, Corr
for i = 1:n_pairs
    r1 = row_idx(i);
    r2 = col_idx(i);
    % Map FC index -> HO label id via roi_ids
    id1 = roi_ids(min(r1, numel(roi_ids)));
    id2 = roi_ids(min(r2, numel(roi_ids)));
    % Label text (fallback if XML missing or index out of range)
    if ~isempty(ho_labels) && id1 >= 1 && id1 <= numel(ho_labels) && ~isempty(ho_labels{id1})
        lab1 = ho_labels{id1};
    else
        lab1 = sprintf('ROI-%d', id1);
    end
    if ~isempty(ho_labels) && id2 >= 1 && id2 <= numel(ho_labels) && ~isempty(ho_labels{id2})
        lab2 = ho_labels{id2};
    else
        lab2 = sprintf('ROI-%d', id2);
    end

    roi_pair_matrix{i,1} = r1;             % index in FC
    roi_pair_matrix{i,2} = lab1;
    roi_pair_matrix{i,3} = r2;             % index in FC
    roi_pair_matrix{i,4} = lab2;
    roi_pair_matrix{i,5} = fc_matrix(r1, r2);
end

roi_pair_table = cell2table(roi_pair_matrix, ...
    'VariableNames', {'ROI_1','Label_1','ROI_2','Label_2','Correlation'});


% -------------------- save it locally --------------------
try
    if opt.SaveTable && ~isempty(opt.SaveDir)
        writetable(roi_pair_table, fullfile(opt.SaveDir, sprintf('%shigh_fc_pairs.csv', prefix)));
    end
catch ME
    warning(ME.identifier, 'Could not save high-FC CSV: %s', ME.message);
end

end
