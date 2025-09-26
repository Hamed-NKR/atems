# TEM Dataset of Wildfire Smoke Particles (Kamloops, BC – July 2021)

This repository contains Transmission Electron Microscopy (TEM) images, processed MATLAB .mat data, and analysis code associated with the study of wildfire smoke particles sampled during the 2021 Southern Interior BC forest fires.  

The dataset complements the ATEMS framework (Automated TEM Segmentation) and uses a forked version maintained here: Hamed-NKR/atems (https://github.com/Hamed-NKR/atems).  

---

## Repository Contents

- images/ : TEM images of biomass smoke particles sampled near Kamloops, BC.
- data/ : MATLAB .mat files containing processed particle data extracted from the TEM images.  
- code/ : MATLAB scripts for segmentation, feature extraction, and visualization. Requires ATEMS.

---

## Installation & Dependencies

1. Install MATLAB (R2021a or later recommended).  
2. Clone the ATEMS repository (or use the fork with wildfire smoke extensions):  

   git clone https://github.com/Hamed-NKR/atems.git

3. Add ATEMS and this dataset repository to your MATLAB path:  

   addpath(genpath('path/to/atems'));
   addpath(genpath('path/to/this_repo'));

---

## Data Description

### Methods

- Particles segmented from TEM images using k-means clustering on pixel intensity values.  
- Boundaries refined and morphological metrics extracted.  
- Pipeline follows Sipkens & Rogak (2021, J. Aerosol Sci.) and Sipkens et al. (2024, JOSS).

### Workflow

The dataset was generated using the following workflow:

1. **Load TEM images** into MATLAB.
2. **Segmentation** using k-means clustering (default) to separate particles from background. Otsu thresholding is also supported but optional.
3. **Binary particle analysis** to extract morphological properties.
4. **Manual classification** of particle morphologies into categories such as tarballs, softballs, soot, ash, minerals, or miscellaneous organics.
5. **Export** of results into MATLAB `.mat` files and summary tables.

### .mat Files
Each .mat file contains:
- pars : table/struct of particle properties  
  - area : projected particle area after segmentation [nm²]  
  - da :  projected area equivalent circular diameter [nm]  
  - ca : circularity (0–1, with 1 = perfect circle)  
  - z_opt : relative optical depth  
  - s_opt : sharpness (boundary acutance)  
  - center_mass [x,y] : centroid coordinates [px,px]  
  - type : classification (tarball, softball, soot, ash, hybrid, etc.)  
  - seg : segmentation masks.  
  - image : original TEM image.
  - binary : binary image after particle segmentation.
and more...

---

## Usage Examples

   \[~, imgs, pixsizes\] = load_imgs
        could also use load_imgs_v3 which can manually extract image scales
   imgs_binary = agg.seg(imgs)  
   pars = analyze_binary_HN(imgs_binary, pixsizes, imgs)
   \[pars, imgs_binary\] = categorize_manu(imgs, pixsizes)  

---

## How to Reproduce Results
1. Use the provided .mat files if you want to skip segmentation.  
2. Or re-run segmentation with provided images and code/ scripts (requires ATEMS).  
3. Compare manual vs. automated classification using boundary sharpness, circularity, and optical depth.
    may also use provided functions such as "spoptics.m" or "categorize_auto.m" to assist with visualization

---

## Citation

If you use this dataset or code, cite:

- Nikookar, H. et al. Image processing methods for morphological characterization of biomass smoke particles, AAAR 2022.  
- Sipkens, T. A. & Rogak, S. N. (2021). Using k-means to identify soot aggregates in TEM images. J. Aerosol Sci. 157, 105807.  
- ATEMS: Sipkens, T. A. et al. (2024). Automated segmentation of TEM images of aerosol particles. JOSS, doi:10.21105/joss.06416.  

---

## License

- Dataset: CC BY 4.0 (or BY-NC 4.0).  
- Code: MIT License (see ATEMS).  

---

## Acknowledgements

We thank Dr. Jason Olfert (UofA), Dr. Allan Bertram (UBC), Dr. Joel Corbin (NRC), and Dr. Arash Naseri (UofA) for their intellectual and logistical support. TEM imaging was performed at UBC Bioimaging Facilities.
