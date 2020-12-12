# ATEMS

**(Matlab *A*nalysis tools for *TEM* images of *S*oot)**

[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)
[![Version](https://img.shields.io/badge/Version-1.0+-blue.svg)]()
[![tsipkens](https://circleci.com/gh/tsipkens/atems.svg?style=shield)]() 

<p align="left">
  <img width="310" src="docs/header.png">
</p>


This codebase contains Matlab code for several methods of characterizing soot aggregates in TEM images. This includes methods for evaluating the aggregate projected area, perimeter, and primary particle diameter. Methods include Otsu thresholding, the pair correlation method (PCM), Hough circle transform (following [Kook](kook)), and tools to aid in manual analysis. This code is designed to replace a [previous, deprecated code](https://github.com/unatriva/UBC-PCM).

**Testing** of this codebase makes use of the `main_*` functions in the upper directory, which are described in [1. MAIN SCRIPTS (main_\*)](#1-main-scripts-main_) below. Specifically, `main_kmeans` and `main_auto` test the fully automated methods, while `main_0` allows for testing of the more manual methods (which require substantial user input). 

The program is primarily composed of two analysis packages, which will be discussed later in the README: 

1. **+agg** - which performs aggregate-level segmentation to output a binary image, and 
2. **+pp** - which determines the primarily particle, often from the binary image generated by the methods in the **agg** package noted above. 

## Table of contents

[Dependencies](#dependencies)

[Getting started](#getting-started): Walking through a sample code

[Components](#1-main-scripts-main_)

- [1. Main scripts (main_\*)](#1-main-scripts-main_): Testing the codebase

- [2. Aggregate segmentation package (+agg)](#2-aggregate-segmentation-package-agg)

- [3. Primary particle analysis package (+pp)](#3-primary-particle-analysis-package-pp)

- [4. Additional tools package (+tools)](#4-additional-tools-package-tools)

[Contributors and acknowledgements](#contributors-and-acknowledgements)

[How to cite](#how-to-cite)

[References](#references)

## Dependencies

This software was tested using Matlab 2020a (though most functions have been validated against older versions) and depends on the following Matlab toolboxes: 

1. the curve fitting toolbox, 
2. the financial toolbox, 
3. the image toolbox, 
4. the optimization toolbox, and 
5. the video and image blockset. 

If not already installed, these packages can be added to Matlab via the instructions available on [Matlab's website](https://www.mathworks.com/help/matlab/matlab_env/get-add-ons.html). 

This program also has an *optional* submodule that adds support for perceptually uniform colormaps compiled from various sources and available at https://github.com/tsipkens/cmap. It can be added via git bash by

```bash
git submodule add -b master https://github.com/tsipkens/cmap
git submodule init
```

The submodule is not necessary for any of the scripts or methods included with this code. 

## Getting started

#### Load images

The first step in the image analysis process is to import images. Images can be handled in one of two ways. The first is as a Matlab structure, with one entry in the structure for each image loaded. To generate this structure by loading the sample images provided with the program, one can use the `tools.load_imgs` function as follows:

```Matlab
Imgs = tools.load_imgs('images');  % load the images
```

The output structure contains the image file name and directory; the image itself, with the footer cropped; and the pixel size, read from the image footer. 

> NOTE: The latter two operations make use of the `tools.detect_footer_scale` function, which attempts to interpret the image footer and/or scale using image analysis tools. The method is known to work with TEM images taken at the University of British Columbia, where it applies optical character recognition to determine the pixel size from the footer text (stored in the `Imgs.pixsize` field) and crops the footer away. Black text on the image may also be readable more generally. Alternatively, the user can also call the `tools.ui_scale_bar` function to use a UI to get the pixel size.

One can also get the image to analyzer by selecting files them in a file explorer by omitting the string input: 

```Matlab
Imgs = tools.load_imgs;  % load the images
```

Proceeding, the image can also be handled is using a cell array of cropped images and pixel sizes. These are secondary outputs to the `tools.load_imgs` function: 

```Matlab
[Imgs, imgs, pixsizes] = tools.load_imgs('images'); % load the images
```

The images and pixel sizes can equivalently be extracted from the `Imgs` structure using:

```Matlab
imgs = {Imgs.cropped}; % copy variables locally
pixsizes = [Imgs.pixsize]; % pixel size for each image
fname = {Imgs.fname};
```

> NOTE: At the moment, this approach will load all of the images into memory. For computers with less memory, this could restrict the number of images that can be analyzed at one time. Batches of 250 images have been run successfully on a computer with 16 GB of RAM (Matlab limits memory usage to below this value). A simple work around is to load the images in groups, rather than all at once, and then add an outer loop to proceed through the different groups. The `tools.load_imgs(fd, n)` function is equipped to handle this, passing the second, optional argument `n` as the integer range of images to load into memory. For example, `n = 1:3` will load the first three images from the set specified by `fd`. 

#### Aggregate-level segmentation

The next step is to evaluate binaries that separate the image into pixels that are part of the background and pixels that are part of aggregates. This is done using the functions in the **agg** package. For example, a *k*-means segmentation specific to this problem can be performed using:

```Matlab
imgs_binary = agg.seg_kmeans(imgs, pixsizes, opts);
    % segment aggregates
```

The result, `imgs_binary`, is either a single binary image (if one image is given as an input) or a cell array of image binaries (if multiple images are given as an input), with `1` if a pixel is considered part of the aggregate and `0` if it is not. Other approaches are also available and are discussed in Section 1 below. 

Having segmented the image, aggregate characteristics can be determined by passing this binary image to an analysis function:

```Matlab
Aggs = agg.analyze_binary(...
    imgs_binary, imgs, pixsizes, fname);
        % determine aggregate properties
```

The output data itself, `Aggs`, is a MATLAB structure with one entry per aggregate, containing properties like the location of the aggregate, perimeter, and projected-area. This data can then be exported to a JSON file using :

```Matlab
tools.write_json(Aggs, fname);
```

or an Excel spreadsheet using: 

```Matlab
tools.write_excel(Aggs, fname);
```

to be analyzed in other software and languages, where `fname` is the filename of the file to write to. The `Aggs` structure can also be visualized as an overlay on top of the TEM image using

```Matlab
figure(1);
tools.imshow_agg(Aggs);
```

The resultant image highlights pixels that are part of the aggregate, and plots a circle that corresponds to the radius of gyration. This output is similar to that shown in images above for the *k*-means and the manual slider methods above. 

#### Determining the primary particle size

Primary particle size information can finally be determined using functions in the **pp** package. The original pair correlation method (PCM), for example, can be applied by using the `Aggs` output from the `agg.analyze_binary` function as

```Matlab
Aggs = pp.pcm(Aggs); % apply PCM
```

The output is an updated `Aggs` structure that also contains an estimate of the primary particle size for each aggregate. Having done this, the primary particle size can be visualized along with the radius of gyration noted above by passing the updated `Aggs` structure to the `tools.plot_aggregates` function:

```Matlab
figure(1);
tools.imshow_agg(Aggs);
```

The inner circle in this plot now indicates the primary particle size from PCM, the larger circle the radius of gyration from *k*-means, and the number indicating the index used by the program to identify each aggregate. Images produced using this type of procedure will feature heavily in the remainder of this README. 

## 1. MAIN SCRIPTS (main_\*)

The main scripts demonstrate further use of the code for specific scenarios and are provided to **test the codebase**. Three such scripts are included: 

**1.** The `main_kmeans` script is designed to specifically investigate the *k*-means approach to aggregate-level segmentation. By default, this is done on the sample images provided in the `images/` folder included with this distribution. This script will also read in some binaries (provided in the `images/slider/` folder) produced by the slider method beforehand for comparison. Finally, the primary particle size is computed for the *k*-means binaries and is plotted for the user. This script supports a European Aerosol Conference submission ([Sipkens et al., 2020][eac20]).

**2.** The `main_auto` script runs through most of the fully automated techniques provided with this program (e.g., *k*-means, Otsu, PCM, etc.), applying them to the images provided in the `images/` folder. The binary images produced by each method are overlaid on the original images for inspection by the user. Outputs are similar to those presented in [Section 2](#2-aggregate-segmentation-package-agg) below. 

**3.** The `main_0` script presents use of the `agg.seg` general segmentation function, described [below](#a-general-segmentation-function-seg), as well as how to read in one's own images. At the end of the script, data is written to an Excel file for examination external to Matlab. Note that this requires that the user's images be formatted appropriately, either with a horizontal footer below the TEM image that has a white background or a cropped TEM image, in which case the user will have to provide pixel size information. Sample images on which the method have been tested are again found in the  `images/` folder. 

## 2. AGGREGATE SEGMENTATION PACKAGE (+agg)

This package contains an expandable library of functions aimed at performing semantic segmentation of the TEM images into aggregate and background regions. Functions are accessed by appending `agg.` to the function name, e.g., `agg.seg_kmeans` to call the *k*-means segmentation procedure. 

We will demonstrate the segmentation functions using a set of size sample images. 

![original](docs/original.png)

These images are included with this distribution in the `images/` folder. These images represent soot collected from a lab-scale flare ([Trivanovic et al., 2020][triv20]) and a diesel engine ([Kheirkhah et al., 2020][kheirkhah20]). 

### 2.1 agg.seg* functions

The main functions implementing aggregate-level semantic segmentation have filenames following `agg.seg_*`. In each case, the output primarily consists of binary images of the same size as the original image but where pixels taken on logical values: `1` for pixels identified as part of the aggregate `0` for pixels identified as part of the background. The functions often take similar inputs, with main argument being `imgs`, which is one of:

1. a single image (after any footer or additional information have been removed), 
2. a cell array of images (again with the footer removed); or
3. an `Imgs` structure, containing the image information as fields. 

Several methods also take `pixsize`, which denotes the size of each pixel in the image. If an `Imgs` structure if provided, this information is expected to be contained in this structure and any `pixsize` input is ignored. Other arguments depend on the function (e.g., optional parameters for the rolling ball transform). 

The set of available methods is summarized below. 

#### A general segmentation function: seg

The `agg.seg` function is a general, multipurpose wrapper function that attempts several methods listed here in sequence, prompting the user after each attempt. Specifically, the method attempts: 

1. The ***k*-means** classifier (following [Sipkens and Rogak][jaskmeans]), prompting the user after segmentation is complete. The user can either (*i*) accept the result as is, (*ii*) reject the output altogether and move on to the next method, (*iii*) choose to remove particles, or (*iv*) add either entirely new particles or add pixels to existing particles (which involves skipping ahead to Step 3, using the current segmentation as a starting point). 

2. The **Otsu** classifier, with the same prompts following segmentation. 

3. Use the GUI-based slider method (the `seg_slider` function described below) to produce a largely-manual segmentation. 

This is repeated until the user has classifid all of the images that were passed to the function. 

#### *k*-means segmentation: seg_kmeans

This function applies a *k*-means segmentation approach following [Sipkens and Rogak][jaskmeans] and using three feature layers, which include: 

FEATURE 1. a *denoised* version of the image, 

FEATURE 2. a measure of *texture* in the image, and 

FEATURE 3. an Otsu-like classified image, with the *threshold adjusted* upwards. 

Compiling these different feature layers results in a three layer image that will be used for segmentation. This is roughly equivalent to segmenting colour images, if each of the layers was assigned a colour. For example, compilation of these feature layers for the sample images results in the following feature layers and compiled RGB image: 

![fcolour](docs/fcolour.png)

Finally, applying Matlab's `imsegkmeans` function, we achieve segmentations as follows: 

![kmeans](docs/kmeans.png)

This is the most robust of the fully automated methods available with this distribution. However, while this will likely be adequate for many users, the technique still occasionally fails, particularly if the function does not adequately remove the background. The method also has some notable limitations when images are (i) *zoomed in* on a single aggregate while (ii) also slightly overexposed. 

#### Otsu thresholding: seg_otsu_rb*

These automated methods apply Otsu thresholding followed by a rolling ball transformation. Two versions of this function are included. 

**1.** The `agg.seg_otsu_rb_orig` function remains more true to the original code of [Dastanpour et al. (2016)][dastanpour2016]. For the sample images, the following segmentations.

![rb_orig](docs/otsu_rb_orig.png)

Aggregates are often broken apart, which may be insufficient in itself. This implementation can be complimented with `agg.seg_slider`, described below, to fill in the gaps between the aggregates and add missing aggregates. 

**2.** Stemming from the deficiencies of the above function, the `agg.seg_otsu_rb` function updates the above implementation by (*i*) not immediately removing boundary aggregates, (*ii*) adding a background subtraction step using the `agg.bg_subtract` function, and (*iii*) adding a bilateral denoising step. This results in the following segmentations.

![otsu_rb](docs/otsu_rb.png)

This latter function generally performs better, though the results still often breaks up aggregates and should likely be compliment with some manual adjustments following initial thresholding. The technique generally underperformed relative to the previously mentioned *k*-means method. 

#### GUI-based slider method: seg_slider

The function `agg.seg_slider` is a largely manual technique. The function enacts a GUI-based method with a slider for adaptive, semi-automatic thresholding of the image (*adaptive* in that small sections of the image can be cropped and assigned individually-selected thresholds). This is done in several steps: 

1. The user is first prompted to **crop** around a single aggregate. This allows the user to zoom in on the image for the remaining steps. 
2. The user uses a lasso tool to draw around the aggregate. The excluded regions are used for **background subtraction** in cropped region of the image. 
3. Gaussian blurring is then performed on the image to reduce the noise in the output binary image. Then, the user is prompted with a **slider** that is used to manually adjust the level of the threshold in the cropped region of the image. Very dark regions denotes sections above the selected threshold. The optimal threshold normally occurs when parts of the surrounding region are also considered above the threshold (i.e., show up as black) but do not connect to the main aggregate (see Step 4 in the figure below depicting the window progression). 
4. The user is prompted to **select** which regions are aggregate, ignoring any white regions that may be above the threshold but are not part of the aggregate. 
5. Finally, the user will be prompted to **check** whether the segmentation was successful by referring to an image with the resultant binary overlaid on the original image.

The progression of windows generally appears as follows. 

![manual_progress](docs/manual_progress.png)

It is worth noting that the mostly manual nature of this approach will resulting in variability and subjectiveness between users. However, the human input often greatly improves the quality of the segmentations and, while more time-intensive, can act as a reference in considering the appropriateness of the other segmentation methods. A sample segmentation follows. 

![slider](docs/slider.png)

Several sub-functions are included within the main file. This is a variant of the method included with the distribution of the PCM code by [Dastanpour et al. (2016)][dastanpour2016]. 

> We note that this code saw important bug updates since the original code by [Dastanpour et al. (2016)][dastanpour2016]. This includes fixing how the original code would repeatedly apply a Gaussian filter every time the user interacted with the slider in the GUI (which may cause some backward compatibility issues), a reduction in the use of global variables, memory savings, and other performance improvements. 

### 2.2 analyze_binary

All of the above methods produce a common output: a binary image. The `agg.analyze_binary` function is now used to convert these binaries to aggregate characteristics, such as area in pixels, radius of gyration, area-equivalent diameter, aspect ratio, among other quantities. The function itself takes a binary image, the original image, and the pixel size as inputs, as follows. 

```Matlab
Aggs = agg.analyze_binary(imgs_binary,imgs,pixsize,fname);
```

The output is a MATLAB structured array, `Aggs`, containing information about the aggregate. The array has one entry for each aggregate found in the image, which is itself defined as any independent groupings of pixels. The `fname` argument is optional and adds this tag to the information in the output `Aggs` structure. 

### 2.3 rolling_ball

Multiple of these methods make use of the *rolling ball transformation*, applied using the `agg.rolling_ball` function included with this package. This transform fills in gaps inside aggregates, acting as a kind of noise filter. This is accomplished by way iterative morphological opening and closing. 

## 3. PRIMARY PARTICLE ANALYSIS PACKAGE (+pp)

The +pp package contains multiple methods for determining the primary particle size of the aggregates of interest. Often this requires a binary mask of the image that can be generated using the +agg package methods. After applying most of the methods, the primary particle size will be stored in (*i*) the `Aggs.dp` field and (*ii*) another `Aggs` field with additional information specifying the method used (e.g., `Aggs.dp_pcm1`  contains a PCM-derived primary particle diameter). For the former, whichever method was last used will overwrite the `Aggs.dp` field, which then acts as a default value that is used by other functions (by the `tools.imshow_agg` function). 

#### PCM

The `pp.pcm` function contains code for the University of British Columbia's pair correlation method (PCM) method. This package contains a significant improvements to code readability, memory use, and length relative to the previous code provided with [Dastanpour et al. (2016)][dastanpour2016]. The underlying method is largely unchanged, correlating the relationship between pixels to the primary particle size for a given aggregate. A single average primary particle diameter is given for each aggregate. [Dastanpour et al.][dastanpour2016] provided two different types of pair correlation function (PCF), corresponding to `Aggs.dp_pcm1`, previously denoted as *simple* and `Aggs.dp_pcm2`, previously denoted as *general*. Testing has generally suggested that the simple method perform better, and this value is assigned to `Aggs.dp` in the output from the PCM method. 

#### EDM-SBS

The Euclidean distance mapping, scale-based analysis (EDM-SBS) of [Bescond et al. (2014)][bescond] is implement in the `pp.edm_sbs` function. This is an adaptation of the original code for use with Matlab and using the binaries above in the place of output from imageJ. As such, some minor differences in output should be expected (which are challenging to compare, as the ImageJ output does not have a direct analog here). The method remains true to how it is described in [Bescond et al.][bescond] and ports some components from the original Scilab code (version 3, available [here](http://www.coria.fr/spip.php?article910)). 

Among the changes to the original EDM-SBS code, this implementation also applies the EDM-SBS method to individual aggregates. While this allows for a better comparison to the other methods here, the method was originally intended to evaluate the primary particle size distribution across a range of aggregates, with uncertaintainty ramifications. However, the linear nature of the curves generated by the method means that the overall EDM-SBS curve is simply approximated by the superposition of all of the aggregates. This is also output by the present code.  

#### Hough transform: kook*

Two `pp.kook*` functions are included with this program, which fit circles to features in the image using the Hough transform and the pre-processing steps described by [Kook et al. (2015)][kook]. 

The first function, `pp.kook`, contains a copy of the code provided by [Kook et al. (2015)][kook], with minor modifications to match the input/output of some of the other packages — namely to take a single image, `img`, and a pixel size, `pixsize` — and to output a `Pp` structure, which contains information for each circle. Note that the original function acts on images without trying to assign primary particles to an aggregate, something resolved in the second function below. This causes some compatibility issues in terms of comparing the output from this function to the other methods contained in this program. 

The `pp.kook2` function contains a modified version of the method proposed by [Kook et al. (2015)][kook] that excludes circles in the background and assigns primary particles to aggregates. This is done rather simply:  by checking if the center of the circles from the preceding procedure lie within the binary for a given aggregate. A sample output is as follows. 

![kook2](docs/kook2.png)

Here, red circles are identified as part of an aggregate, while black circles are excluded in the output. Note that it is apparent from some of the images that gradients in the background influenced the results. Though, it can also be noted that the circles within the aggregates remain reasonable (to the extent that the overall method is reasonable). 

#### Manual sizing

The `pp.manual` function can be used to manual draw circles around the soot primary particles. The code was developed at the University of British Columbia and represents a heavily modified version of the code associated with [Dastanpour and Rogak (2014)][dastanpour2014]. The current method uses two lines overlaid on each primary particle to select the length and width of the particle. This is converted to various quantities, such as the center of each primary particle and the overall mean primary particle diameter. 


## 4. ADDITIONAL TOOLS PACKAGE (+tools)

This package contains a series of functions that help in visualizing or analyzing the aggregates that transcend multiple methods or functions. We discuss a few examples here and refer the reader to individual functions for more information. 

#### Functions to show images (tools.imshow*)

These functions aid in visualizing the resultant images. For example, 

```Matlab
tools.imshow_binary(img, img_binary);
```

will plot the image, given in `img`, and overlay labels for the aggregates in a corresponding binary image, given in `img_binary`.  Appending an `opts` structure to the function allows for the control of the appearance of the overlay, including the label alpha and colour. For example, 

```Matlab
opts.cmap = [0.92,0.16,0.49];
tools.imshow_binary(img, img_binary, opts);
```

will plot the overlays in a red, while 

```Matlab
opts.cmap = [0.99,0.86,0.37];
tools.imshow_binary(img, img_binary, opts);
```

will plot the overlays in a yellow. 

The related `tools.imshow_agg` function will plot the binaries, as above, and add aggregate-level information to the plot, including: (*i*) the aggregate number; (*ii*) the radius of gyration about the center of the aggregate; and (*iii*) the average primary particle diameter, if this information is available (the function will use the `Aggs.dp` field, which will contain primary particle information from the most recently applied method). 

#### Functions to write data to files (tools.write*)

This set of functions writes data to files, with the precise format depending on the function. 

The `write_excel` and `write_json` functions write aggregate-level information to Excel and JSON formats respectively. The precise output will depend on what information is contained in the `Aggs` structure given to the method. For example, `Aggs` structures with primary particle information will have that information output to the Excel or JSON files. 

The `write_images` function, in contrast, take a cell of images (or a single image) and writes them to a series of files, specified by the `fnames` argument, in the folder specified by the `folder` argument. For example, this can be used to write the binary images presented at the beginning of this README to a temporary folder using

```Matlab
tools.imwrite(imgs_binary, fname, 'temp');
```

This function is also useful when paired with the `tools.imshow_agg` function to write images with binary overlays and radius of gyration information (such as those presented throughout this README). For a series of images that have been processed to produce an `Aggs` structure, this can be accomplished using: 

```Matlab
imgs_agg{length(imgs)} = []; % initialize cell

% Loop through images and get overlay for each image. 
% The second argument of tools.imshow_agg gets the 
% colordata from the frame in the overlay image. 
for ii=1:length(imgs)
    [~, imgs_agg{ii}] = tools.imshow_agg(Aggs, ii, 1, opts);
end

% Write the overlay images to a temporary folder
tools.write_images(imgs_agg, fname, 'temp');
```

As the image size will depend on the physical size of the figure on the screen, it may be useful to first run

```Matlab
f1 = figure(1); f1.WindowState = 'maximized';
```

to maximize the figure window, to generate higher quality output. 

#### Functions to visualize post-processed results (tools.viz*)

This set of functions is intended to generate formatted plots of post-processed results. For example, `tools.viz_darho` generates a scatter plot of the area-equivalent diameter versus effective density estimated for each aggregate. 

--------------------------------------------------------------------------

#### License

This software is relaesed under an MIT license (see the corresponding license file for details).

#### Contributors and acknowledgements

This code was primarily compiled by Timothy A. Sipkens while at the University of British Columbia (UBC), who can be contacted at [tsipkens@mail.ubc.ca](mailto:tsipkens@mail.ubc.ca). 

Pieces of this code were adapted from various sources and features snippets written by several individuals at UBC, including Ramin Dastanpour, [Una Trivanovic](https://github.com/unatriva), Alberto Baldelli, Yiling Kang, Yeshun (Samuel) Ma, and Steven Rogak, among others.

This program contains very significantly modified versions of the code distributed with [Dastanpour et al. (2016)][dastanpour2016]. The most recent version of the Dastanpour et al. code prior to this overhaul is available at https://github.com/unatriva/UBC-PCM (which itself presents a minor update from the original). That code forms the basis for some of the methods underlying the manual processing and the PCM method used in this code, as noted in the README above. However, significant optimizations have improved code legibility, performance, and maintainability (e.g., the code no longer uses global variables). 

Also included with this program is the Matlab code of [Kook et al. (2015)][kook], modified to accommodate the expected inputs and outputs common to the other functions.

This code also contain an adaptation of EDM-SBS method of [Bescond et al. (2014)][bescond]. We thank the authors, in particular Jerome Yon, for their help in understanding their original [Scilab code and ImageJ plugin](http://www.coria.fr/spip.php?article910). Modifications to allow the method to work directly on binary images (rather than a custom output from ImageJ) and to integrate the method into the Matlab environment may present some minor compatibility issues, but allows use of the aggregate segmentation methods given in the **agg** package. 

Finally, the progress bar in the function `tools.textbar`, which is used to indicate progress on some of the primary particle sizing techniques, is a modified version of that written by [Samuel Grauer](https://github.com/sgrauer). 

#### How to cite

On use of this code or the *k*-means segmentation procedure described [above](#k-means-segmentation-seg_kmeans), please cite: 

> [Sipkens, T. A., Rogak, S. N. (Accepted). Technical note: Using k-means to identify soot aggregates in transmission electron microscopy images. *Journal of Aerosol Science*, 105699.][jaskmeans]

Users of the pair correlation method (PCM), the Euclidean distance mapping-surface based scale (EBD-SBS), and Hough transform (following Kook et al.) codes should acknowledge the corresponding studies under the acknowledgements above. Please also consider citing this repository directly. 

#### References

[Bescond, A., Yon, J., Ouf, F. X., Ferry, D., Delhaye, D., Gaffié, D., Coppalle, A. & Rozé, C. (2014). Automated determination of aggregate primary particle size distribution by TEM image analysis: application to soot. Aerosol Science and Technology, 48(8), 831-841.][bescond]

[Dastanpour, R., Boone, J. M., & Rogak, S. N. (2016). Automated primary particle sizing of nanoparticle aggregates by TEM image analysis. Powder Technology, 295, 218-224.][dastanpour2016]

[Dastanpour, R., & Rogak, S. N. (2014). Observations of a correlation between primary particle and aggregate size for soot particles. Aerosol Science and Technology, 48(10), 1043-1049.][dastanpour2014]

[Kheirkhah, P., Baldelli, A., Kirchen, P. & Rogak, S., (2020). Development and validation of a multi-angle light scattering method for fast engine soot mass and size measurements. Aerosol Science and Technology, 54(9), 1083-1101.][kheirkhah20]

[Kook, S., Zhang, R., Chan, Q. N., Aizawa, T., Kondo, K., Pickett, L. M., Cenker, E., Bruneaux, G., Andersson, O., Pagels, J., & Nordin, E. Z. (2016). Automated detection of primary particles from transmission electron microscope (TEM) images of soot aggregates in diesel engine environments. *SAE International Journal of Engines*, *9*(1), 279-296.][kook]

[Sipkens, T. A., Zhou., L., Rogak, S. N. (2020). Aggregate-level segmentation of soot TEM images by unsupervised machine learning. *European Aerosol Conference*, Aachen, Germany.][eac20]

[Sipkens, T. A., Rogak, S. N. (Accepted). Technical note: Using k-means to identify soot aggregates in transmission electron microscopy images. *Journal of Aerosol Science*, 105699.][jaskmeans]

[Trivanovic, U., Sipkens, T. A., Kazemimanesh, M., Baldelli, A., Jefferson, A. M., Conrad, B. M., Johnson, M. R., Corbin, J. C., Olfert, J. S., & Rogak, S. N. (2020). Morphology and size of soot from gas flares as a function of fuel and water addition. *Fuel*, *279*, 118478.][triv20]

[kook]: https://doi.org/10.4271/2015-01-1991
[dastanpour2016]: https://doi.org/10.1016/j.powtec.2016.03.027
[dastanpour2014]: https://doi.org/10.1080/02786826.2014.955565
[bescond]: https://doi.org/10.1080/02786826.2014.932896
[triv20]: https://doi.org/10.1016/j.fuel.2020.118478
[eac20]: https://doi.org/10.13140/RG.2.2.14433.12648
[jaskmeans]: https://doi.org/10.1016/j.jaerosci.2020.105699
[kheirkhah20]: https://doi.org/10.1080/02786826.2020.1758623
