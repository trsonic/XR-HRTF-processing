# XR-based HRTF Measurement Postprocessing
This set of MATLAB functions is a part of the [XR-based HRTF Measurement System](https://trsonic.github.io/XR-HRTFs/). It serves as a postprocessing workflow for binaural sine sweep signals recorded using the [HRTF Measurement Control App](https://github.com/trsonic/XR-HRTF-capture).

## Getting Started
* Requirements
    * MATLAB (tested using R2022a)
    * [SOFA MATLAB Toolbox](https://github.com/sofacoustics/SOFAtoolbox)
    * [Ambisonic Decoder Toolbox](https://bitbucket.org/ambidecodertoolbox/adt/src/master/) - optionally, comment out `saveAsAmbix` function call if you don't need Ambisonic decoders to be calculated.
* Usage
    * Once the binaural recordings have been captured, edit the path to the subject dir in `IRprocessing.m` script. Running the script should execute all necessary processing.

## Output data
* A set of various plots (located in the subject folder under `/figures/`)
* SOFA files:
    * RAW HRIRs at measured directions
    * Diffuse-field Equalized HRIRs at measured directions
    * RAW HRIRs at interpolated directions
    * Diffuse-field Equalized HRIRs at interpolated directions
* Measured HRIRs in Wave format
* Ambix Ambisonic decoder config files
