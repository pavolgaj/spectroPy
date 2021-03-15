# spectroPy

![](https://img.shields.io/github/languages/top/pavolgaj/spectropy.svg?style=flat)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/b82d1a45153c41169fcff464b15d4924)](https://app.codacy.com/gh/pavolgaj/spectroPy?utm_source=github.com&utm_medium=referral&utm_content=pavolgaj/spectroPy&utm_campaign=Badge_Grade_Settings)
![GitHub all releases](https://img.shields.io/github/downloads/pavolgaj/spectropy/total?label=GitHub&nbsp;downloads)

 Calibration and manipulation with grating spectra in Python. Focused mainly on spectra obtained using StarAnalyser SA-100 and SA-200.
 
 __Requirements:__
 * python
 * numpy
 * scipy
 * matplotlib
 * astropy
 
 
 __Recommended workflow:__
 
 1. Take images of spectrum.
 2. Apply dark-frame and flat-field corrections, register and stack images etc. (optional).
 3. Extract L channel (optional).
 4. Extract line profile of image with 0th order a 1st order spectrum - use e.g. SAO DS9 or MaximDL.
 5. Start spectroPy: `python spectroPy.py your_profile`.
 6. Run calibration of spectrum:
   * Use 1-point calibration if you know a spectral resolution of used equipment.
   * Use 2-points calibration otherwise.
   * If your image was color and there is three channels in your profile, you can run calibration for only one channel or for all channels separately (slower but more precise).
 7. Select region of profile with 0th order spectrum - use left mouse button to mark left border of this region and right button for right border.
 8. Check the fitted gauss profile. If it looks bad, try to select slightly different region of spectra (bigger or smaller).
 9. If you choosed 1-point calibration, add a spectral resolution in Angstrom per pixel.
 10. If you choosed 2-points calibration, select area with known spectral line (similarly to steps 7 and 8). And write the wavelength of this spectral line (in Angstrom).
 11. Now, you can crop your spectrum, plot it and save it to a CSV file, text file or FITS file.
 
