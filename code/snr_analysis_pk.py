from astropull_pk import ImagePull
import numpy as np


def arrays(pos, line):
    """
    I think this function can be separated as it parses your catalogue input. We want to try to have astropull
    independent of catalogue formats, etc., so best to get the various source properties and pass them to
    ImagePull.
    """
    list = []
    num = line.shape[0]
    for n in range(num):
        list.append(line[n][pos])
    list = np.asarray(list)

    return list


def parse_cat_file(filename, src_num):
    """
    Wrap the fetching of columns from the file in a function

    Input
    =======================
    filename (str):
        catalogue filename

    src_num (int):
        the number of the source in the catalogue
    """
    # Loads csv file
    line = np.genfromtxt(open(filename, "r"), names=True, delimiter=',', dtype=None)

    # Arrays Function
    # MCSNR
    MCSNR = arrays(0, line)
    MCSNR = MCSNR[src_num]
    MCSNR = MCSNR.decode("utf-8")

    # RA
    RA = arrays(1, line)
    RA = RA[src_num]

    # DE
    DE = arrays(2, line)
    DE = DE[src_num]

    # Radius
    Rad = arrays(8, line)
    Rad = Rad[src_num]

    # kT
    kT = arrays(12, line)
    kT = kT[src_num]

    # VShock
    VShock = arrays(16, line)
    VShock = VShock[src_num]

    # Age
    Age = arrays(18, line)
    Age = Age[src_num]

    # LX
    LX = arrays(4, line)
    LX = LX[src_num]

    # LIR
    LIR = arrays(7, line)
    LIR = LIR[src_num]

    return MCSNR, RA, DE, Rad, kT, VShock, Age, LX, LIR


# get SNR properties
SNR = 34 # position in list
MCSNR, RA, DE, Rad, kT, VShock, Age, LX, LIR = parse_cat_file('SNR_list.csv', SNR)
instrument_list = ['IRAC', 'MIPS']
sensor = 'IRAC'
b = 8.0
my_test = ImagePull(MCSNR, RA, DE, Rad, sensor, b)
my_test.run()
#for SNR in range(len(np.genfromtxt(open('SNR_list.csv', "r"), names=True, delimiter=',', dtype=None))):
    

#    MCSNR, RA, DE, Rad, kT, VShock, Age, LX, LIR = parse_cat_file('SNR_list.csv', SNR)


#    # instrument things
#    instrument_list = ['IRAC', 'MIPS']

    # Select Instrument and Band
#    for sensor in instrument_list:
#        if sensor == 'IRAC':
#            band_list = [3.6, 4.5, 5.8, 8.0]
#            for b in band_list:
#                # Run image pull as object
#                my_test = ImagePull(MCSNR, RA, DE, Rad, sensor, b)
#                my_test.run()

#        if sensor == 'MIPS':
#            band_list = [24., 70., 160.]
#            for b in band_list:
#                # Run image pull as object
#                my_test = ImagePull(MCSNR, RA, DE, Rad, sensor, b)
#                my_test.run()