class ImagePull:

    def arrays(self):
        # make usable for all arrays
        MCSNR = []
        num = line.shape[0]
        for n in range(num):
            name = line[n][0]
            MCSNR.append(name)
        MCSNR = np.asarray(MCSNR)
        return MCSNR

    def ra2dec(self):
        # debug for function
        x = [Ra.decode('UTF-8') for Ra in RA]
        y = [Dec.decode('UTF-8') for Dec in DE]

        # Creates Ra and Dec formats from data
        im_ra = []
        for n in range(num):
            things = x[n].split()
            Ra = things[0] + 'h' + things[1] + 'm' + things[2] + 's'
            im_ra.append(Ra)

        im_dec = []
        for n in range(num):
            things = y[n].split()
            Dec = things[0] + 'd' + things[1] + 'm' + things[2] + 's'
            im_dec.append(Dec)

        return im_ra, im_dec

    def query(self):
        my_query = sha.query(coord=coord.SkyCoord(ra=test_src_coord[0], dec=test_src_coord[1],
                                                  unit=(u.degree, u.degree)), size=0.001)

        sensor_l = len(sensor) + 1
        sens = sensor.rjust(sensor_l)
        name = str(sens) + ' ' + str(wavelength) + 'um'
        if sensor == 'MIPS':
            if wavelength >= 100:
                name_l = len(name)
            else:
                name_l = len(name) + 1
            name = name.ljust(name_l)
        print(name)
        mask = (my_query['wavelength'] == name)
        if sensor == 'IRAC':
            filesize_mask = (my_query['filesize'] > 2.e8)  # keep only files >200MB (i.e., SAGE images)
            combined_mask = np.logical_and(mask, filesize_mask)
        else:
            combined_mask = mask

        tbl = my_query[combined_mask]

        nimages = len(tbl)
        distance_check = np.zeros(nimages)

        if sensor == 'IRAC':
            for n, row in enumerate(tbl):
                im_ra = row['ra']
                im_dec = row['dec']

                # determine the angular separation between the image centre and SNR
                snr = SkyCoord(test_src_coord[0] * u.degree, test_src_coord[1] * u.degree, frame='icrs')
                im_centre = SkyCoord(im_ra * u.degree, im_dec * u.degree, frame='icrs')
                sep = snr.separation(im_centre)
                distance_check[n] = sep.value

            closest = np.argmin(distance_check)
            my_final_row = tbl[closest]
            print(my_final_row['modedisplayname', 'wavelength', 'filetype', 'filesize', 'ptcomment'])
            print("\n")
            print("Getting: {}".format(my_final_row['accessUrl']))
            url = my_final_row['accessUrl'].strip()
            sha.save_file(url, out_dir='test_files' + '_' + str(wavelength) + '/')
        elif sensor == 'MIPS':
            for i in tbl['accessUrl']:
                print("Getting: {}".format(i))
                url = i.strip()
                sha.save_file(url, out_dir='test_files' + '_' + str(wavelength) + '/')

        return

    def cuts(self):

        my_image_files = glob.glob(os.path.join('test_files_' + str(wavelength), '*.fits'))
        n_images = len(my_image_files)
        counter = 0
        filename = sens[1:] + '_' + str(wavelength) + '.fits'
        array = []
        # cycle through and cut images
        for i in my_image_files:
            with fits.open(i) as hdulist:
                header = hdulist[0].header
                if sensor == 'MIPS':
                    if header['OBJECT'] != 'LMC' or header['OBSRVR'] != 'Margaret Meixner':
                        if os.path.exists(os.path.join(my_image_files[counter])):
                            os.remove(os.path.join(my_image_files[counter]))

        for i in range(len(my_image_files)):
            if not os.path.exists(os.path.join(MCSNR[SNR])):
                os.makedirs(os.path.join(MCSNR[SNR]))
            if os.path.exists(os.path.join(my_image_files[i])):
                copyfile(os.path.join(my_image_files[i]), os.path.join((MCSNR[SNR].decode("utf-8")) + '/' + filename))

        rmtree(os.path.join('test_files' + '_' + str(wavelength) + '/'))

        return

    def plot_image(data, wcs, coords=None):
        """
        convience function to plot data on a wcs projection

        Input:
        data -    the data array
        wcs -     the WCS object
        coords -  the coordinates of the object for plotting (1x2 array, optional)
        """
        # set up the plot
        fig, axs = plt.subplots(1, 1, figsize=(8, 8), subplot_kw={'projection': wcs})

        # set up limits of the colour scale for plotting
        vmin = 1e-2

        # plot
        axs.imshow(data, cmap='jet', interpolation='nearest', origin='lower', norm=LogNorm(vmin=vmin))
        axs.scatter(coords[0], coords[1], transform=axs.get_transform('fk5'), s=300, lw=3,
                    edgecolor='white', facecolor='none')
        axs.set_facecolor('black')
        axs.coords.grid(True, color='white', ls='dotted')
        axs.coords[0].set_axislabel('Right Ascension (J2000)')
        axs.coords[1].set_axislabel('Declination (J2000)')

        # display the plot
        plt.tight_layout()
        plt.show()

    def im_stack(self):

        image_concat = []
        for image in my_image_files:
            image_concat.append(fits.getdata(image))

        final_image = np.zeros(shape=image_concat[0].shape)

        for image in image_concat:
            final_image += image

        # todo impliment drizzle
        # outfile = 'stacked_' + MCSNR[SNR] + '.fits'

        # hdu = fits.PrimaryHDU(final_image)
        # hdu.writeto(outfile, clobber=True)

        plt.imshow(final_image, cmap='jet')

        return

    def ap_phot(self):

        # cycle through images
        # def ap_phot
        fluxes = []
        fluxesav = []
        im_name = []
        for i in my_image_files:
            # load the file data, header, and wcs
            with fits.open(i) as hdulist:
                my_hdu = hdulist[0]
                my_hdu.data = np.nan_to_num(my_hdu.data)
                r = Rad[SNR] * u.arcsec

                position = SkyCoord(test_src_coord[0] * u.degree, test_src_coord[1] * u.degree, frame='icrs')
                apertures = SkyCircularAperture(position, r=Rad[SNR] * u.arcsec)

                # print(apertures)
                phot_table = aperture_photometry(my_hdu, apertures)

                # todo debug
                # onverted_aperture_sum = (phot_table['aperture_sum'] * arc_pix)
                fluxes.append(phot_table['aperture_sum'])
                # im_name.append(my_image_files(i))

                x = Angle(Rad[SNR], u.arcsec)
                y = 2 * x.degree

                t1_pos = SkyCoord((test_src_coord[0] + y) * u.degree, (test_src_coord[1] + y) * u.degree, frame='icrs')
                t2_pos = SkyCoord((test_src_coord[0] + y) * u.degree, (test_src_coord[1] - y) * u.degree, frame='icrs')
                t3_pos = SkyCoord((test_src_coord[0] - y) * u.degree, (test_src_coord[1] + y) * u.degree, frame='icrs')
                t4_pos = SkyCoord((test_src_coord[0] - y) * u.degree, (test_src_coord[1] - y) * u.degree, frame='icrs')

                ap1 = SkyCircularAperture(t1_pos, r=Rad[SNR] * u.arcsec)
                ap2 = SkyCircularAperture(t2_pos, r=Rad[SNR] * u.arcsec)
                ap3 = SkyCircularAperture(t3_pos, r=Rad[SNR] * u.arcsec)
                ap4 = SkyCircularAperture(t4_pos, r=Rad[SNR] * u.arcsec)

                phot_table1 = aperture_photometry(my_hdu, ap1)
                fluxesav.append(phot_table1['aperture_sum'])
                phot_table2 = aperture_photometry(my_hdu, ap2)
                fluxesav.append(phot_table2['aperture_sum'])
                phot_table3 = aperture_photometry(my_hdu, ap3)
                fluxesav.append(phot_table3['aperture_sum'])
                phot_table4 = aperture_photometry(my_hdu, ap4)
                fluxesav.append(phot_table4['aperture_sum'])
                average = np.mean(fluxesav)
                flux = fluxes - average
                print(flux)

        return
