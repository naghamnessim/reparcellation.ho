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

2. **Run the Pipeline as illustrated in the video demo**
