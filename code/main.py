#main
from ImagePull_temp import ImagePull
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
MCSNR = ImagePull.arrays(0)
#RA
RA = ImagePull.arrays([1])
#DE
DE = ImagePull.arrays([2])
#Rad
Rad = ImagePull.arrays([8])
#kT
kT = ImagePull.arrays([12])
#VShock
VShock = ImagePull.arrays([16])
#Age
Age = ImagePull.arrays([18])
#LX
LX = ImagePull.arrays([4])
#LIR
LIR = ImagePull.arrays([7])

# Eq2Deg Function
test_src_coord = ImagePull.eq2dec(RA(SNR), DEC(SNR))

# Query Function
ImagePull.query(test_src_coord, sensor, wavelength)

# Cuts Function
ImagePull.cuts()

# Plot Image function
ImagePull.plot_image(data, wcs, coords=test_src_coord)

# Ap Phot Function
flux = ImagePull.ap_phot()
