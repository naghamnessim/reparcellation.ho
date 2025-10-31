# Re-parcellation of ds000224 fMRI Dataset Using the Harvard–Oxford Atlas

This repository contains MATLAB code, scripts, and documentation for re-parcellating
the **OpenNeuro ds000224** fMRI dataset according to the **Harvard–Oxford (HO) atlas**.
The pipeline includes preprocessing, registration to T1 and MNI space,
atlas-based projection, and ROI-level functional connectivity analysis.

---

## Project Structure
**reparcellation.ho/**
1. code/
- preprocessing/
- projection/
- analysis/
- MainPipeline.mlx/
2. docs/
- DATASET_LINKS.md
- VIDEO_LINK.md
- PRESENTATION/
3. README.md

---

## Pipeline Overview

The complete workflow performs:

1. **Preprocessing**  
   - Slice timing correction and motion correction  
   - Skull stripping for both T1 and EPI volumes  
   - Coregistration of functional to anatomical space  

2. **Projection and Normalization**  
   - T1 normalization to MNI space (`normalize_T1_to_mni.m`)  
   - Warping of functional images to MNI (`warp_epi_to_mni.m`)  
   - Projection of the Harvard–Oxford atlas to native EPI space  
   - Visualization of ROIs overlaid on the EPI volume  

3. **Analysis**  
   - Extraction of ROI time series  
   - Computation of functional connectivity (correlation matrices)  
   - Reporting of high-correlation ROI pairs  

---

## Requirements

- **MATLAB** R2022a or later  
- **SPM12** installed and added to your MATLAB path  
- Sufficient disk space for the ds000224 dataset

---

## Usage

1. **Download the Dataset**

   Download **ds000224 (version 00002)** from OpenNeuro:  
   [https://openneuro.org/datasets/ds000224/versions/00002](https://openneuro.org/datasets/ds000224/versions/00002)

   Place it inside your project directory.

2. **Run the Pipeline as illustrated in the video demo [open here](docs/VIDEOS_LINK.md)**

---

## Demo Results on sub-MSC01 (ses-func02) using the memorywords task

<!-- 1) Underlays -->
<p align="center">
  <img src="docs/Example on Results (sub-MSC01)/HO_EPI Underlay.png" width="48%">
  <img src="docs/Example on Results (sub-MSC01)/HO_T1 Underlay.png" width="48%">
</p>

*Figure 1. Harvard–Oxford atlas underlay views in native EPI and T1 spaces.*

<!-- 2) Preprocessing check (mean EPI before/after skull strip) -->
<p align="center">
  <img src="docs/Example on Results (sub-MSC01)/step5_mean_epi_pre_strip.png" width="48%">
  <img src="docs/Example on Results (sub-MSC01)/step5_mean_epi_post_strip.png" width="48%">
</p>

*Figure 2. Mean EPI before vs. after skull stripping.*

<!-- 3) Atlas overlay on mean EPI -->
<p align="center">
  <img src="docs/Example on Results (sub-MSC01)/step5_overlay.png" width="70%">
</p>

*Figure 3. Harvard–Oxford atlas projected to native EPI (overlay on mean EPI).*

<!-- 4) Motion parameters (if you want to show QC) -->
<!-- Note: filename contains an apostrophe; see note below. -->
<p align="center">
  <img src="docs/Example on Results (sub-MSC01)/step1_motion_params'.png" width="70%">
</p>

*Figure 4. Motion parameters (QC).*

<!-- 5) Functional connectivity -->
<p align="center">
  <img src="docs/Example on Results (sub-MSC01)/step9_fc_matrix.png" width="70%">
</p>

*Figure 5. Parcel-wise functional connectivity matrix.*

<!-- 6) ROI time series overview -->
<p align="center">
  <img src="docs/Example on Results (sub-MSC01)/step9_roi_timeseries_overview.png" width="70%">
</p>

*Figure 6. ROI time-series overview (z-scored).*

<!-- 7) High-FC pairs table link (optional) -->
**Data table:** [Top high-correlation ROI pairs (`step10_high_fc_pairs.csv`)](docs/Example%20on%20Results%20%28sub-MSC01%29/step10_high_fc_pairs.csv)

---
