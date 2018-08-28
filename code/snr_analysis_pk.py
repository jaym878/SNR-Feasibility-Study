from astropull_pk import ImagePull
from astropy import units as u
import numpy as np
import os

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from astropy.table import Table

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


# setup master output array
col_names = ('name', 'radius', 'kT', 'v_s', 'age', 'Lx',
             '3.6_flux', '3.6_arc_flux', '3.6_Jy', '3.6_arc_Jy',
             '4.5_flux', '4.5_arc_flux', '4.5_Jy', '4.5_arc_Jy',
             '5.8_flux', '5.8_arc_flux', '5.8_Jy', '5.8_arc_Jy',
             '8.0_flux', '8.0_arc_flux', '8.0_Jy', '8.0_arc_Jy',
             '24_flux', '24_arc_flux', '24_Jy', '24_arc_Jy',
             '60_flux', '60_arc_flux', '60_Jy', '60_arc_Jy',
             '120_flux', '120_arc_flux', '120_Jy', '120_arc_Jy')

col_dtype = ('S10', 'f4', 'f4', 'f4', 'f4', 'f4',
             'f4', 'f4', 'f4', 'f4',
             'f4', 'f4', 'f4', 'f4',
             'f4', 'f4', 'f4', 'f4',
             'f4', 'f4', 'f4', 'f4',
             'f4', 'f4', 'f4', 'f4',
             'f4', 'f4', 'f4', 'f4',
             'f4', 'f4', 'f4', 'f4')

master_output = Table(None, names=col_names, masked=True, dtype=col_dtype)

# determine number of SNRs in catalogue
n_snrs = len(np.genfromtxt(open('SNR_list.csv', "r"), names=True, delimiter=',', dtype=None))

# get SNR properties using ImagePull
for SNR in range(n_snrs):

    MCSNR, RA, DE, Rad, kT, VShock, Age, LX, LIR = parse_cat_file('SNR_list.csv', SNR)

    print("\nRunning {}\n".format(str(MCSNR)))

    # write the known parameters to master output
    snr_row = [MCSNR, Rad, kT, VShock, Age, LX, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0,0, 0, 0, 0,
               0, 0, 0, 0,0, 0, 0, 0,0, 0, 0, 0]

    # instrument things
    instrument_list = ['IRAC', 'MIPS']
    arc_flux_plot = []
    arc_flux_plot_jy = []

    # Select Instrument and Band
    for sensor in instrument_list:
        flux = []
        arc_flux = []
        flux_jy = []
        arc_flux_jy = []

        if sensor == 'IRAC':
            band_list = [3.6, 4.5, 5.8, 8.0]
            for b in band_list:
                try:
                    # Run image pull as object
                    my_test = ImagePull(MCSNR, RA, DE, Rad, sensor, b)
                    f_jy, f_erg, farc_jy, farc_erg = my_test.run()

                    flux.append(f_erg)
                    arc_flux.append(farc_erg)
                    arc_flux_plot.append(farc_erg)
                    flux_jy.append(f_jy)
                    arc_flux_jy.append(farc_jy)
                    arc_flux_plot_jy.append(farc_jy)

                    if b == 3.6:
                        snr_row[6] = f_erg
                        snr_row[7] = farc_erg
                        snr_row[8] = f_jy
                        snr_row[9] = farc_jy
                    if b == 4.5:
                        snr_row[10] = f_erg
                        snr_row[11] = farc_erg
                        snr_row[12] = f_jy
                        snr_row[13] = farc_jy
                    if b == 5.8:
                        snr_row[14] = f_erg
                        snr_row[15] = farc_erg
                        snr_row[16] = f_jy
                        snr_row[17] = farc_jy
                    if b == 8.0:
                        snr_row[18] = f_erg
                        snr_row[19] = farc_erg
                        snr_row[20] = f_jy
                        snr_row[21] = farc_jy

                except Exception as e:
                    print('{}-{}um failed'.format(sensor, b))
                    print('  %s: %s' % (e.__class__.__name__, str(e)))

                    try:
                        os.mkdir(str(MCSNR))
                    except:
                        pass

                    flux.append(np.nan)
                    arc_flux.append(np.nan)
                    arc_flux_plot.append(np.nan)
                    flux_jy.append(np.nan)
                    arc_flux_jy.append(np.nan)
                    arc_flux_plot_jy.append(np.nan)

                    if b == 3.6:
                        snr_row[6] = np.nan
                        snr_row[7] = np.nan
                        snr_row[8] = np.nan
                        snr_row[9] = np.nan
                    if b == 4.5:
                        snr_row[10] = np.nan
                        snr_row[11] = np.nan
                        snr_row[12] = np.nan
                        snr_row[13] = np.nan
                    if b == 5.8:
                        snr_row[14] = np.nan
                        snr_row[15] = np.nan
                        snr_row[16] = np.nan
                        snr_row[17] = np.nan
                    if b == 8.0:
                        snr_row[18] = np.nan
                        snr_row[19] = np.nan
                        snr_row[20] = np.nan
                        snr_row[21] = np.nan



            band_list_arr = np.asarray(band_list)
            flux_arr = np.asarray(flux)
            flux_arc_arr = np.asarray(arc_flux)
            flux_jy_arr = np.asarray(flux_jy)
            flux_arc_jy_arr = np.asarray(arc_flux_jy)

            flux_out = np.column_stack((band_list_arr, flux_arr))
            flux_arc_out = np.column_stack((band_list_arr, flux_arc_arr))

            flux_jy_out = np.column_stack((band_list_arr, flux_jy_arr))
            flux_arc_jy_out = np.column_stack((band_list_arr, flux_arc_jy_arr))

            out_name = str(MCSNR) + '/fluxes_' + str(MCSNR) + '_' + str(sensor) + '.txt'
            np.savetxt(out_name, flux_out, delimiter=',')

            out_name = str(MCSNR) + '/arc_fluxes_' + str(MCSNR) + '_' + str(sensor) + '.txt'
            np.savetxt(out_name, flux_arc_out, delimiter=',')

            out_name = str(MCSNR) + '/fluxes_Jy_' + str(MCSNR) + '_' + str(sensor) + '.txt'
            np.savetxt(out_name, flux_jy_out, delimiter=',')

            out_name = str(MCSNR) + '/arc_fluxes_Jy_' + str(MCSNR) + '_' + str(sensor) + '.txt'
            np.savetxt(out_name, flux_arc_jy_out, delimiter=',')


        if sensor == 'MIPS':
            band_list = [24, 70, 160]
            for b in band_list:
                try:
                    # Run image pull as object
                    my_test = ImagePull(MCSNR, RA, DE, Rad, sensor, b)
                    f_jy, f_erg, farc_jy, farc_erg = my_test.run()

                    flux.append(f_erg)
                    arc_flux.append(f_erg)
                    arc_flux_plot.append(farc_erg)
                    flux_jy.append(f_jy)
                    arc_flux_jy.append(farc_jy)
                    arc_flux_plot_jy.append(farc_jy)

                    if b == 24:
                        snr_row[22] = f_erg
                        snr_row[23] = farc_erg
                        snr_row[24] = f_jy
                        snr_row[25] = farc_jy
                    if b == 70:
                        snr_row[26] = f_erg
                        snr_row[27] = farc_erg
                        snr_row[28] = f_jy
                        snr_row[29] = farc_jy
                    if b == 160:
                        snr_row[30] = f_erg
                        snr_row[31] = farc_erg
                        snr_row[32] = f_jy
                        snr_row[33] = farc_jy

                except Exception as e:
                    print('{}-{}um failed'.format(sensor, b))
                    print('  %s: %s' % (e.__class__.__name__, str(e)))

                    try:
                        os.mkdir(str(MCSNR))
                    except:
                        pass

                    flux.append(np.nan)
                    arc_flux.append(np.nan)
                    arc_flux_plot.append(np.nan)
                    flux_jy.append(np.nan)
                    arc_flux_jy.append(np.nan)
                    arc_flux_plot_jy.append(np.nan)

                    if b == 24:
                        snr_row[22] = np.nan
                        snr_row[23] = np.nan
                        snr_row[24] = np.nan
                        snr_row[25] = np.nan
                    if b == 70:
                        snr_row[26] = np.nan
                        snr_row[27] = np.nan
                        snr_row[28] = np.nan
                        snr_row[29] = np.nan
                    if b == 160:
                        snr_row[30] = np.nan
                        snr_row[31] = np.nan
                        snr_row[32] = np.nan
                        snr_row[33] = np.nan



            band_list_arr = np.asarray(band_list)
            flux_arr = np.asarray(flux)
            flux_arc_arr = np.asarray(arc_flux)
            flux_jy_arr = np.asarray(flux_jy)
            flux_arc_jy_arr = np.asarray(arc_flux_jy)

            flux_out = np.column_stack((band_list_arr, flux_arr))
            flux_arc_out = np.column_stack((band_list_arr, flux_arc_arr))

            flux_jy_out = np.column_stack((band_list_arr, flux_jy_arr))
            flux_arc_jy_out = np.column_stack((band_list_arr, flux_arc_jy_arr))

            out_name = str(MCSNR) + '/fluxes_' + str(MCSNR) + '_' + str(sensor) + '.txt'
            np.savetxt(out_name, flux_out, delimiter=',')

            out_name = str(MCSNR) + '/arc_fluxes_' + str(MCSNR) + '_' + str(sensor) + '.txt'
            np.savetxt(out_name, flux_arc_out, delimiter=',')

            out_name = str(MCSNR) + '/fluxes_Jy_' + str(MCSNR) + '_' + str(sensor) + '.txt'
            np.savetxt(out_name, flux_jy_out, delimiter=',')

            out_name = str(MCSNR) + '/arc_fluxes_Jy_' + str(MCSNR) + '_' + str(sensor) + '.txt'
            np.savetxt(out_name, flux_arc_jy_out, delimiter=',')

    # add everything to the master output table
    master_output.add_row(snr_row)

    # plot in SED erg/cm2/s
    fig, axs = plt.subplots(1, 1, figsize=(5, 4))

    band_list = [3.6, 4.5, 5.8, 8.0, 24, 70, 160]

    axs.plot(band_list, arc_flux_plot, 'r-', marker='o', markersize=3, label="avg-flux/fov")
    name = str(MCSNR) + '-Flux vs Wavelength'
    axs.set_title(name)
    axs.set_xlabel('Wavelength (um)')
    axs.set_ylabel(r'Flux (erg s$^{-1}$ cm$^{-2}$)')
    axs.set_xscale('log')
    plt.tight_layout()

    plot_name = os.path.join(str(MCSNR) + '/flux_plot_' + str(MCSNR) + '.pdf')

    try:
        os.remove(plot_name)
    except:
        pass

    # plot SED in Jy
    fig.savefig(plot_name, dpi=200)

    fig, axs = plt.subplots(1, 1, figsize=(5, 4))

    band_list = [3.6, 4.5, 5.8, 8.0, 24, 70, 160]

    arc_flux_plot_jy_arr = np.asarray(arc_flux_plot_jy)
    axs.plot(band_list, arc_flux_plot_jy_arr * 1.e3, 'r-', marker='o', markersize=3, label="Flux/MIRI FOV")
    name = str(MCSNR) + '-Flux vs Wavelength'
    axs.set_title(name)
    axs.set_xlabel('Wavelength (um)')
    axs.set_ylabel('Flux density (mJy)')
    axs.set_xscale('log')
    plt.tight_layout()

    plot_name = os.path.join(str(MCSNR) + '/flux_Jy_plot_' + str(MCSNR) + '.pdf')

    try:
        os.remove(plot_name)
    except:
        pass

    fig.savefig(plot_name, dpi=200)


# save the master output to file
master_name = 'SNR-catalogue_IR_properties.txt'

try:
    os.remove(master_name)
except:
    pass

master_output.write(master_name, format='ascii')
