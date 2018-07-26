from astropull import ImagePull
import numpy as np

# Loads csv file
line = np.genfromtxt(open('/home/murphyj/Desktop/Coding/SNR_list.csv', "r"), names=True, delimiter=',', dtype=None)

SNR = 34 # position in list
instrument_list = ['IRAC', 'MIPS']
band_list = [3.6, 4.5, 5.8, 8.0, 24, 70, 160]
# Select Instrument and Band
sensor = instrument_list[0]
wavelength = band_list[1]

# Arrays Function
#MCSNR
MCSNR = ImagePull.arrays(0, line)
#RA
RA = ImagePull.arrays(1, line)
#DE
DE = ImagePull.arrays(2, line)
#Rad
Rad = ImagePull.arrays(8, line)
#kT
kT = ImagePull.arrays(12, line)
#VShock
VShock = ImagePull.arrays(16, line)
#Age
Age = ImagePull.arrays(18, line)
#LX
LX = ImagePull.arrays(4, line)
#LIR
LIR = ImagePull.arrays(7, line)

# Eq2Deg Function
test_src_coord = ImagePull.eq2deg(RA[SNR], DE[SNR])

# Query Function
ImagePull.query(test_src_coord, sensor, wavelength)

# Cuts Function
ImagePull.cuts(wavelength, sensor, MCSNR[34])

# Plot Image function
# cycle through images
my_image_files = ImagePull.show_image(MCSNR[SNR], test_src_coord)

# Ap Phot Function
flux = ImagePull.ap_phot(my_image_files, Rad[SNR], test_src_coord, wavelength)
print(flux)
