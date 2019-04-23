ImageJ Scripts: Use to pre-process raw .tif files obtained from imaging rounds
  1. Make maximum intensity projection
  2. Crop
  3. Split into individual channels
  
STORM Analysis: check https://github.com/ZhuangLab/storm-analysis

MATLAB: To process data obtained from STORM analysis
  Takes .HDF5 files
  Uses oligodT-A488 to mask and define cell boundary, DAPI to define nuclear boundary
  .xslx files are used as placeholders for the MATLAB script to write the results into
