function overlay_ROIs_volshow(func_or_struct_path, atlas_path, volume_index)
% Inputs:
%   func_or_struct_path (char) — path to an underlay NIfTI volume readable by SPM (3D or 4D EPI/T1)
%   atlas_path          (char) — path to a labeled ROI atlas (integer labels, 0=background) in the same grid
%   volume_index        (int, optional) — which 3D frame to show if the underlay is 4D (e.g. 60)
%
% Outputs:
%   (none) — launches MATLAB’s volumeViewer with the underlay + ROI labels overlaid
%
% Purpose:
%   Display a 3D anatomical (or a single 3D frame from a 4D fMRI) with a parcel label map overlaid.
%   Uses robust window/level via percentiles and gentle gamma to enhance visibility, and crops to the
%   ROI bounding box for a focused view and faster rendering. The atlas is shown with nearest-neighbor
%   categorical overlay inside volumeViewer.
%
% Author: Nagham Nessim
% University of Geneva, 2025

    if nargin < 3, volume_index = 60; end

    % ---------- Load underlay (3D, or a chosen frame from 4D) ----------
    % spm_vol reads the header(s); spm_read_vols reads voxel data in double.
    V = spm_vol(func_or_struct_path);
    if numel(V) > 1
        vol = spm_read_vols(V(volume_index));  % Extract requested 3D frame from 4D
    else
        vol = spm_read_vols(V);                % already 3D frame 
    end

    % ---------- Load atlas (labeled integer volume; background = 0) ----------
    atlas_img = spm_read_vols(spm_vol(atlas_path));

    % ---------- Basic grid match check and force NaNs to background for windowing ----------
    if ~isequal(size(vol), size(atlas_img))
        warning('overlay_ROIs_volshow:GridMismatch', ...
            'Underlay size %s vs atlas size %s. Ensure the atlas was resliced onto the underlay grid correctly.', ...
            mat2str(size(vol)), mat2str(size(atlas_img)));                           % Mismatch will cause overlay to not align voxelwise.
    end

    vol(isnan(vol)) = min(vol(~isnan(vol)));

    % ---------- Compute ROI bounding box to crop both volumes ----------
    roi_mask = atlas_img > 0;
    if ~any(roi_mask(:))
        warning('Atlas has no nonzero labels; opening full volume without overlay.');
        % Window/level + gamma for the full volume
        p     = prctile(vol(:), [2 98]);             % define a window
        vol_w = mat2gray(vol, [p(1) p(2)]);          % rescale to [0,1]
        gamma = 0.8;                                 % Gamma<1 brightens midtones; helps soft-tissue contrast
        try
            vol_w = imadjustn(vol_w, [], [], gamma); % Image Processing Toolbox is needed for this step
        catch
            vol_w = rescale(vol_w .^ gamma, 0, 1);   % fallback if imadjustn is unavailable
        end
        vol_uint8 = uint8(255 * vol_w);
        volumeViewer(vol_uint8);                     % launch volumeViewer
        return
    end
    
    % Crop to labeled bounding box
    [x, y, z] = ind2sub(size(roi_mask), find(roi_mask));
    x1 = min(x); x2 = max(x);
    y1 = min(y); y2 = max(y);
    z1 = min(z); z2 = max(z);

    vol_crop = vol(x1:x2, y1:y2, z1:z2);
    roi_crop = atlas_img(x1:x2, y1:y2, z1:z2);

    % ---------- Window/level + gentle gamma on the cropped underlay ----------
    % - Windowing via robust percentiles (2–98) improves contrast while ignoring outliers.

    p = prctile(vol_crop(:), [2 98]);
    vol_w = mat2gray(vol_crop, [p(1) p(2)]);
    gamma = 0.8;
    try
        vol_w = imadjustn(vol_w, [], [], gamma);
    catch
        vol_w = rescale(vol_w .^ gamma, 0, 1);
    end

    vol_uint8 = uint8(255 * vol_w);

    % ---------- Display with volumeViewer ----------
    try
        volumeViewer(vol_uint8, roi_crop);
    catch ME
        error('overlay_ROIs_volshow:ViewerError', ...
              'volumeViewer failed: %s. Ensure you have the Image Processing Toolbox (volumeViewer).', ME.message);
    end
end

