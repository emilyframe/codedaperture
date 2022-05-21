# Introduction
`codedaperture` is a python library for reconstructing gamma-ray images from coded aperture imager data. This library contains scripts for building the coded aperture system response ('sysresponse.py'), performing MLEM image reconstruction ('mlem.py'), and providing a graphical image of the reconstructed image data ('showimage.py').

# Python Dependencies
- numpy
- tables
- sys
- time
- matplotlib.pyplot (if using the 'showimage.py' script)
- matplotlib.animation (if using the 'showimage.py' script)

# Scripts ('codedaperture/scripts')
- `sysresponse.py`
  - Builds the system response and saves the system response to `codedaperture/sysresp/sysresp<pose #>.h5`
  - The system response is stored as a 2D matrix of i x j elements, where i is the number source voxels in the image space and j is the number of detector pixels.
  - The system response matrix can be produced at various projection angles (i.e. detector poses) for tomographic image reconstruction. If multiple detector poses are given in the `poses` array, a system response matrix will be saved to a separate file for each pose and named accordingly.
  - Note: this script reads in data from `codedaperture/efficiency.h5` and `codedaperture/farfieldmap.h5`. You will need to download `farfieldmap.h5` from: https://drive.google.com/file/d/1CcjEMqjV6u81_Q0JffcZDwmgxi3MNiyR/view?usp=sharing and store this file in `codedaperture/`

- `mlem.py`
  - Performs MLEM image reconstruction using the system response data in `codedaperture/sysresp/sysresp<pose #>.h5` and the detector data in `codedaperture/sysresp/data/spiral<pose #>.h5`.
  - This script allows for tomographic image reconstruction, i.e. combining system response matrices and detector data at multiple detector poses.
    - Note: be sure to match system response matrix with corresponding detector data
  - The reconstructed data is saved to `codedaperture/images/image.h5`.

- `showimage.py`
  - Converts the image data in `codedaperture/images/image.h5` to a 3D volume and animates 2D slices of this volume over various depths (perpendicular to the imager).

# Instructions for Running the Scripts
1. In the `codedaperture` directory, run the command `mkdir sysresp` and `mkdir images`.
1. Build system response matrix: `python scripts/sysresponse.py`.
  - Note: change data pathway at top of `codedaperture/scripts/sysresponse.py` to the location at which you stored this repo.
  - Note: it takes roughly 3 minutes to build one system response.
2. Perform MLEM image reconstruction: `python scripts/mlem.py`
  - Note: change data pathway at top of `codedaperture/scripts/mlem.py` to the location at which you stored this repo.
3. Plot reconstructed data: `python scripts/showimage.py`
  - Note: change data pathway at top of `codedaperture/scripts/showimage.py` to the location at which you stored this repo.
