from astropull_pk import ImagePull
import numpy as np

if __name__ == "__main__":
    MCSNR = sys.argv[0]
    RA = sys.argv[1]
    DE = sys.argv[2]
    Rad = sys.argv[3]

#    # instrument things
    instrument_list = ['IRAC', 'MIPS']

    # Select Instrument and Band
    for sensor in instrument_list:
        if sensor == 'IRAC':
            band_list = [3.6, 4.5, 5.8, 8.0]
            for b in band_list:
                # Run image pull as object
                my_test = ImagePull(MCSNR, RA, DE, Rad, sensor, b)
                my_test.run()

        if sensor == 'MIPS':
            band_list = [24., 70., 160.]
            for b in band_list:
                # Run image pull as object
                my_test = ImagePull(MCSNR, RA, DE, Rad, sensor, b)
                my_test.run()
