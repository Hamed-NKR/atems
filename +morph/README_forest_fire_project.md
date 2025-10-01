# TEM Dataset of Wildfire Smoke Particles (Kamloops, BC – July 2021)

**Authors**  
Hamed Nikookar<sup>a</sup>, Timothy A. Sipkens<sup>b</sup>, Steven N. Rogak<sup>a</sup>  

**Affiliations**  
<sup>a</sup> Department of Mechanical Engineering, University of British Columbia,  
6250 Applied Science Lane, Vancouver, V6T 1Z4, British Columbia, Canada  

<sup>b</sup> Metrology Research Centre, National Research Council Canada,  
1200 Montreal Road, Ottawa, K1A 0R6, Ontario, Canada  

**Contact**  
[hmdnkr@student.ubc.ca]

This repository contains Transmission Electron Microscopy (TEM) images, processed MATLAB `.mat` data, and analysis code associated with the study of wildfire smoke particles sampled during the 2021 Southern Interior BC forest fires.  

The dataset complements the `ATEMS` framework (Automated TEM Segmentation; `https://github.com/tsipkens/atems`) and uses a forked version maintained here: `Hamed-NKR/atems` (`https://github.com/Hamed-NKR/atems`).  

---

## Repository Contents

- `images/` : TEM images of biomass smoke particles sampled in Kamloops, BC, collected during daytime and nighttime.
- `data/` : MATLAB `.mat` files containing processed particle data extracted from the TEM images.  
- `codes/` : MATLAB scripts for segmentation, feature extraction, and visualization. Requires `ATEMS`.

---

## Installation & Dependencies

1. Install MATLAB (R2021a or later recommended).  
2. Clone the `ATEMS` repository (or use the fork with wildfire smoke extensions):  

   `git clone https://github.com/Hamed-NKR/atems.git;`

3. Add `ATEMS` and this dataset repository to your MATLAB path:  

   `addpath(genpath('path/to/atems'));`

   `addpath(genpath('path/to/this_repo'));`

---

## Data Description

### Methods

- Particles segmented from TEM images using k-means clustering on pixel intensity values.  
- Morphological metrics extracted, particularly: projected area, circularity, optical depth, and sharpness.
- Visually decided particle type lables.
- Pipeline follows Sipkens & Rogak (2021, J. Aerosol Sci.), Nikookar et al. (2022, AAAR), and Sipkens et al. (2024, JOSS).

### Workflow

The dataset was generated using the following workflow:

1. **Load TEM images** into MATLAB.
2. **Segmentation** using k-means clustering (default) to separate particles from background. Otsu thresholding is also supported but optional.
3. **Binary particle analysis** to extract morphological properties.
4. **Manual classification** of particle morphologies into categories such as tarballs, softballs (rounded organics), soot, ash-like, minerals, hybrids (type mixtures) or miscellaneous.
5. **Export** of results into MATLAB `.mat` files and summary tables.

### .mat Files
Each .mat file contains:
- `pars` : table/struct of particle properties  
  - `area` : projected particle area after segmentation [nm²]  
  - `da` :  projected area equivalent circular diameter [nm]  
  - `ca` : circularity (0-1, with 1 = perfect circle)  
  - `z_opt` : relative optical depth (0-1, with 1 = perfectly black)
  - `s_opt` : sharpness (boundary acutance, higher implies sharper boundaries)  
  - `center_mass` [x,y] : centroid coordinates [px,px]  
  - `type` : classification (tarball, softball, soot, ash, hybrid, etc.)  
  - `image` : original TEM image.
  - `binary` : segmentation masks.
and more...

---

## Usage Examples

1. `[~, imgs, pixsizes] = load_imgs`

   May, alternatively, use `load_imgs_v3.m` which can manually extract image scales.

2. `imgs_binary = agg.seg(imgs)`  
3. `pars = analyze_binary_HN(imgs_binary, pixsizes, imgs)`
4. `[pars, imgs_binary] = categorize_manu(imgs, pixsizes)`  

---

## How to Reproduce Results
1. Use the provided .mat files if you want to skip segmentation.  
2. Or re-run segmentation with provided images and code/ scripts (requires `ATEMS`).  
3. Compare manual vs. automated classification using boundary sharpness, circularity, and optical depth.

   May use provided functions such as `spoptics.m` or `categorize_auto.m` to assist with visualization.

---

## Citation

If you use this dataset or code, cite:

- Nikookar, H. et al. Image processing methods for morphological characterization of biomass smoke particles, AAAR 2022, [doi:10.5281/zenodo.15588129](https://zenodo.org/records/15588129).  
- ATEMS: Sipkens, T. A. et al. (2024). Automated segmentation of TEM images of aerosol particles. JOSS 9(99), 6416, [doi:10.21105/joss.06416](https://doi.org/10.21105/joss.06416).  
- Sipkens, T. A. & Rogak, S. N. (2021). Using k-means to identify soot aggregates in TEM images. J. Aerosol Sci. 152, 105699, [doi:10.1016/j.jaerosci.2020.105699](https://doi.org/10.1016/j.jaerosci.2020.105699).  

---

## License

- Dataset: CC BY 4.0.  
- Code: GNU General Public License v3.0.  

---

## Acknowledgements

We thank Dr. Jason Olfert (UofA), Dr. Allan Bertram (UBC), Dr. Joel Corbin (NRC), and Dr. Arash Naseri (UofA) for their intellectual and logistical support. Samples were collected by Mang Guan (UBC) and imaged/analyzed by the first author. Thermophoretic sampler for particle collection was developed by Steven Zimmerman and Muhammad M. H. Zareer (UBC). The data processing was conducted under supervision of Dr. Steven Rogak and Dr. Timothy Sipkens. TEM imaging was performed at UBC Bioimaging Facilities.
