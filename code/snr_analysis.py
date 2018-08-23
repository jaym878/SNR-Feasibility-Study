from astropull import ImagePull
from astropy import units as u
import numpy as np
import os

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

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
#SNR = 4 # position in list
for SNR in range(len(np.genfromtxt(open('/home/murphyj/Desktop/Coding/SNR_list.csv', "r"), names=True, delimiter=',', dtype=None))):

    MCSNR, RA, DE, Rad, kT, VShock, Age, LX, LIR = parse_cat_file('/home/murphyj/Desktop/Coding/SNR_list.csv', SNR)


    # instrument things
    instrument_list = ['IRAC', 'MIPS']
    h_flux_val = []
    l_flux_val = []
    # Select Instrument and Band
    for sensor in instrument_list:
        flux = []
        
        if sensor == 'IRAC':
            band_list = [3.6, 4.5, 5.8, 8.0]
            for b in band_list:
                try:
                    # Run image pull as object
                    my_test = ImagePull(MCSNR, RA, DE, Rad, sensor, b)
                    ret = my_test.run()
                    FluxIR = ret[0]
                    x = ret[1] * (1 * u.MJy)
                    flux_high = x[0].value
                    x = ret[2] * (1 * u.MJy)
                    flux_low = x[0].value

                    flux.append(FluxIR)
                    h_flux_val.append(flux_high)
                    l_flux_val.append(flux_low)

                except Exception:
                    print(Warning: missing data)
                    FluxIR = 0
                    flux_high = 0
                    flux_low = 0
                    flux.append(FluxIR)
                    h_flux_val.append(flux_high)
                    l_flux_val.append(flux_low)

            with open('flux_template_irac.txt', 'w') as filehandle:  
                filehandle.writelines("%s\n" % place for place in band_list)
            with open(os.path.join(str(MCSNR) + '/fluxes_' + str(MCSNR) + '_' + str(sensor) + '.txt'), 'w') as filehandle:  
                filehandle.writelines("%s\n" % place for place in flux)
            with open(os.path.join(str(MCSNR) + '/high_fluxes_' + str(MCSNR) + '_' + str(sensor) + '.txt'), 'w') as filehandle:  
                filehandle.writelines("%s\n" % place for place in h_flux_val)
            with open(os.path.join(str(MCSNR) + '/low_fluxes_' + str(MCSNR) + '_' + str(sensor) + '.txt'), 'w') as filehandle:  
                filehandle.writelines("%s\n" % place for place in l_flux_val)

        if sensor == 'MIPS':
            band_list = [24, 70, 160]
            for b in band_list:
                try:
                    # Run image pull as object
                    my_test = ImagePull(MCSNR, RA, DE, Rad, sensor, b)
                    ret = my_test.run()
                    FluxIR = ret[0]
                    x = ret[1] * (1 * u.MJy)
                    flux_high = x[0].value
                    x = ret[2] * (1 * u.MJy)
                    flux_low = x[0].value

                    flux.append(FluxIR)
                    h_flux_val.append(flux_high)
                    l_flux_val.append(flux_low)

                except Exception:
                    print(Warning: missing data)
                    FluxIR = 0
                    flux_high = 0
                    flux_low = 0
                    flux.append(FluxIR)
                    h_flux_val.append(flux_high)
                    l_flux_val.append(flux_low)

            with open('flux_template_irac.txt', 'w') as filehandle:  
                filehandle.writelines("%s\n" % place for place in band_list)
            with open(os.path.join(str(MCSNR) + '/fluxes_' + str(MCSNR) + '_' + str(sensor) + '.txt'), 'w') as filehandle:  
                filehandle.writelines("%s\n" % place for place in flux)
            with open(os.path.join(str(MCSNR) + '/high_fluxes_' + str(MCSNR) + '_' + str(sensor) + '.txt'), 'w') as filehandle:  
                filehandle.writelines("%s\n" % place for place in h_flux_val)
            with open(os.path.join(str(MCSNR) + '/low_fluxes_' + str(MCSNR) + '_' + str(sensor) + '.txt'), 'w') as filehandle:  
                filehandle.writelines("%s\n" % place for place in l_flux_val)

    band_list = [3.6, 4.5, 5.8, 8.0, 24, 70, 160]
    plt.plot(band_list, h_flux_val, 'r', label="High Flux")
    plt.plot(band_list, l_flux_val, 'b', label="Low Flux")
    name = str(MCSNR) + '-Flux vs Wavelength'
    plt.title(name)
    plt.xlabel('Wavelength (um)')
    plt.ylabel('Flux (erg/s/cm^2)')

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    plot_name = os.path.join(str(MCSNR) + '/flux_plot_' + str(MCSNR) + '.pdf')

    plt.savefig(plot_name, dpi=200)