# carboseg

This folder contains Python code to implement a convolutional neural network (CNN) for segmentation as described by [Sipkens et al.][ptech.cnn] Details and code for the training of the network are available in a parallel repository at https://github.com/maxfrei750/CarbonBlackSegmentation, with primary contributions by Max Frei ([@maxfrei750](https://github.com/maxfrei750)). The implementation here makes use of the ONNX file output (to be downloaded [here](https://uni-duisburg-essen.sciebo.de/s/J7bS47nZadg4bBH/download)) from that procedure and employs the Python ONNX runtime for execution. Use of this function requires the necessary Python environment as a pre-requisite. 

### **Setup: Creating the necessary Python environment**

If [conda](https://conda.io/en/latest/miniconda.html) is installed, one can use the following procedure to create the necessary Python environment: 

**1.** First, open the conda command line. 

**2.** Change into the `carboseg/` directory in a cloned copy of the current repository, using:

```shell
cd [path]/carboseg
```

where `[path]` is the path the parent MATLAB repository. 

**3.** Create a new Python environment using the `environment.yml` file included in that directory, using:

```shell
conda env create --file environment.yml
```

This will create the **carboseg** environment and the corresponding Python executable with the necessary dependencies. 

Alternatively, one can set up an environment to take advantage of a GPU. The example here applied to CUDA-enabled scenerios. In this case, one can create a **carboseg-gpu** environmnet using:

```shell
conda env create --file environment-gpu.yml
```

One must then point to the appropriate alternative Python interpretter. Now, one can apply the CNN, either (*i*) directly through MATLAB or (*ii*) by using MATLAB in conjunction with a Python IDE. 

### **Segmentation directly using MATLAB**

To start, the current distribution includes a function to aid in loading the Python environment:

```Matlab
py_exec = path_to_exe;
tools.load_python;
```

where `path_to_exe` is the path to the Python executable for the **carboseg** environment. The user should edit the script to point to an appropriate Python executable generated by the previous procedure  For Windows, this code also adds the necessary folders and a reference to the local [carboseg](https://github.com/tsipkens/atems/tree/master/carboseg) folder. 

Finally, the user can use the `agg.seg_carboseg(...)` function in a similar fashion to the other segmentation approaches, noting the limitations at the beginning of this section. Optimally, the function should run with the standard inputs and outputs common to the other classifiers discussed here. 

For an implementation of this procedure, see the `main_carboseg` script in the upper directory of this repository.

> NOTE: There are known issues with the MATLAB calls to Python, including freezing indefinitely in the `seg_carboseg(...)` function call and Python environment errors. Potential workaround for the freezing/hanging is to add a debug point in the `seg_carbonseg(...)` function near the top of the file or to restart MATLAB and run the scripts in a new session. In these instances, it may be better to save pre-processed images (with footer cropped), run the code in Python, and reload the images into MATLAB. See the next subsection for this option. 

> NOTE: If CUDA is available and the batch of images is reasonably large, it may be faster to save the images and run the classification on a GPU in Python directly (again, see the next subsection for this option). 

### **Segmentation using Python (with read/write to MATLAB)**

Alternatively, one can save the images, load them in a Python function directly, save the results, and reload the images in MATLAB. This can be broken into three steps. First, load the images, as before, and save the cropped images to a folder, `fd_in`, for reading in Python,  

```Matlab
% Load images (here using the 'images' folder).
[Imgs, imgs, pixsizes] = tools.load_imgs('images');
fnames = {Imgs.fname};  % also extract file names

agg.seg_ext(imgs, fnames, fd_in);  % save the images to `fd` folder
```

Second, in a Python IDE, one can now run the `main.py` script in the [carboseg](https://github.com/tsipkens/atems/tree/master/carboseg) folder, editing the script to point to the appropriate folder where MATLAB saved the cropped images from the previous step. Finally, one can return to MATLAB and read in processed images using the `agg.seg_cnn_pt2(...)` method, 

```Matlab
imgs_binary = agg.seg_cnn_pt2(fnames, fd_out, pixsizes)
```

where `fd_out` is the new folder containing the classified binaries from Python. The `imgs_binary` contains a cell of binary images, matching the other classifiers here. Note that the image file name should be made consistent between calls, such that the images can be appropriately matched when transitioning back and forth between MATLAB and Python. One can now proceed with post-processing analogous to the other classifiers. 

For an implementation of this procedure, see the `main_carboseg_ext` script in the upper directory of this repository.

#### Reference

[Sipkens, T.A., Frei, M., Baldelli, A., Kirchen, P., Kruis, F. E., & Rogak, S. N. (2021) Characterizing soot in TEM images using a convolutional neural network. Powder Technology.][ptech.cnn]

[ptech.cnn]: https://doi.org/10.1016/j.powtec.2021.04.026
