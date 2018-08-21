# SNR-Feasibility-Study
### Goal
The SNR Feasibility study was designed to look at supernova remnants (SNRs) in the Large Magellanic Cloud (LMC). The package does this by accessing the Spitzer heritage archive and downloading images in all wavelengths utilised by the IRAC and MIPS instruments. The package analyses these SNRs by referencing uploaded data by the user containing the co-ordinates of the central point and the radius of the SNR. The package then applies cuts and aperture photometry to obtain the average flux of the SNR. 
### How to use
The package requires inputs of the name of the target, equatorial coordinates such as 06 32 56, -05 56 24, the radius in parsecs, the sensor being used and the waveband for testing. This can be used in the form 
```
ImagePull(name, ra, dec, radius, sensor, wavelength)
```
This then returns a 1x3 array of fluxes. The 0th position contains the average flux across the area of the SNR, the 1st position contains the low average flux across the circumference of the SNR and the 2nd position contains the high average flux across the circumference of the SNR. These fluxes are returned with the units erg/s/cm^2.
### Required file shape
Currently the package requires the data to be uploaded by reading in a .csv file where the required data is in the format
```
NAME      | RA      | DEC       | RADIUS
----------------------------------------
J0448-6700| 04 48 22| -66 59 52 | 26.67
```
