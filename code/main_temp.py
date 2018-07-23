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
# I want to pass the "0" to the function in ImagePull_temp.py
MCSNR = ImagePull.arrays(0)
