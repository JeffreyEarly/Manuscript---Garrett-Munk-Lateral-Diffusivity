garrett-munk-lateral-diffusivity/data
===========
This folder contains files that have been extracted from the original model runs.

The mooring time series .mat files were extracted using,
    GLOceanKit/Matlab/WintersModel/Examples/SaveMooringTimeSeries.m

The wave-vortex decomposition can be done using,
    GLOceanKit/Matlab/WintersModel/Examples/SaveWaveVortexDecomposition.m
This creates a very large NetCDF file that is basically the same size as the original model output (it's an equivalent representation). The first time step takes a long time, but it should go more quickly after that.

The analysis of the decorrelation time series is done with,
    GLOceanKit/Matlab/WintersModel/Examples/SaveWaveVortexDecorrelationTimeMatrix.m

