# utilities
#main
import os, glob
from pprint import pprint
from shutil import copyfile, rmtree
from photutils import SkyCircularAperture
from photutils import SkyCircularAnnulus
from photutils import SkyRectangularAperture
from photutils import aperture_photometry
from photutils import datasets

# standard packages
import numpy as np
import math
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

# astropy and astroquery
from astroquery import sha
from astropy import coordinates as coord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy.nddata import Cutout2D


class ImagePull:
    """
    ImagePull references the Spitzer Heritage Archive using supplied coordinates to pull images 
    from both the IRAC and MIPS sensor wavebands. The class then applies cuts and aperture photometry to derive fluxes of targets.
    Write what this class does in the doc string here. Include whatever input it needs and what it outputs

    Input
    ==========================
    name (str):
        name of the SNR

    ra (float):
        ra of the SNR

    dec (float):
        dec of the SNR

    radius (float):
        radius of the SNR in units of arcseconds

    sensor (str):
        Spitzer detector (options are IRAC and MIPS)

    wavelength (float):
        Spitzer filter required (options are: 3.6, 4.5, 5.8, 8.0, 24., 70., 160.)


    Output
    ==========================
    cut images:
        saved pdf images of target in selected wavebands
    
    FluxIR (float):
        flux returned in units of erg/s/cm^2
    """
    def __init__(self, name, ra, dec, radius, sensor, wavelength):
        # load input parameters to attributes
        self.name = name
        self.ra = ra
        self.dec = dec
        self.radius = radius
        self.sensor = sensor
        self.wavelength = wavelength

        self.ra_deg = None
        self.dec_deg = None

        self.flux_Jy = None
        self.flux_erg = None


    def eq2deg(self, RA, DE):
        """
        eq2deg converts a string of eq coordinates to a decimal degree equivelant

        Input
        ==========================
        RA (str):
            ra in units of h/m/s

        DE (str):
            dec in units of d/m/s

        Output
        ==========================
        im_ra_deg (float):
            ra in units of decimal degrees

        im_dec_deg (float):
            dec in units of decimal degrees
        """
        # debug for function
        x = RA.decode('UTF-8')
        y = DE.decode('UTF-8')

        # Creates Ra and Dec formats from data
        things = x.split()
        Ra = things[0] + 'h' + things[1] + 'm' + things[2] + 's'
        im_ra = Ra

        things = y.split()
        Dec = things[0] + 'd' + things[1] + 'm' + things[2] + 's'
        im_dec = Dec

        a = Angle(im_ra, u.hour)
        im_ra_deg = a.degree

        b = Angle(im_dec)
        im_dec_deg = b.degree

        self.ra_deg = im_ra_deg
        self.dec_deg = im_dec_deg

        return im_ra_deg, im_dec_deg

    def query(self, test_src_coord, sensor, wavelength):
        """
        eq2deg converts a string of eq coordinates to a decimal degree equivelant
        
        Input
        ==========================
        test_src_coord (float):
            coordinates of SNR in decimal degrees
            
        sensor (str):
            sensor name
            
        wavelength (str):
            waveband for sensor
            
        """
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
            #print(name)
        mask = (my_query['wavelength'] == name)
        if sensor == 'IRAC':
            filesize_mask = (my_query['filesize'] > 2.e8)
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
                snr = SkyCoord(test_src_coord[0] * u.degree, test_src_coord[1] * u.degree, frame='icrs')
                im_centre = SkyCoord(im_ra * u.degree, im_dec * u.degree, frame='icrs')
                sep = snr.separation(im_centre)
                distance_check[n] = sep.value

            closest = np.argmin(distance_check)
            my_final_row = tbl[closest]
            url = my_final_row['accessUrl'].strip()
            sha.save_file(url, out_dir=str(self.name) + '_' + str(wavelength) + '/')

        elif sensor == 'MIPS':
            for i in tbl['accessUrl']:
                print("Getting: {}".format(i))
                url = i.strip()
                sha.save_file(url, out_dir=str(self.name) + '_' + str(wavelength)+ '/')

    def fill_sources(self, image, save_plot=True):
        """
        detect sources, mask and fill the gaps using photutils

        Input
        =========================

        image (2D np array)
            image to process

        save_plot (boolean)
            if True, saves a diagnostic plot of pre and post source removal

        Returns
        =========================

        filled_image (2D np array)
            image with sources removed and filled
        """
        from photutils.background import Background2D
        from astropy.stats import gaussian_fwhm_to_sigma
        from photutils import detect_sources
        from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans, convolve

        # get the background level
        bkg = Background2D(image, 30)
        threshold = bkg.background + (5.0 * bkg.background_rms)

        # create the size of the sources
        sigma = 2.0 * gaussian_fwhm_to_sigma  # FWHM = 2.
        kernel = Gaussian2DKernel(sigma, x_size=2, y_size=2)
        kernel.normalize()

        # search for the sources and create a mask
        seg_image = detect_sources(image, threshold, npixels=5, filter_kernel=kernel)
        source_mask = np.clip(seg_image, 0, 1)

        # fill the masked source regions
        masked_image = image.copy()
        masked_image[source_mask > 0] = np.nan

        kernel = Gaussian2DKernel(x_stddev=5)
        filled_image = interpolate_replace_nans(masked_image, kernel)

        # smooth out hard edges around filled sources
        kernel = Gaussian2DKernel(x_stddev=3)
        filled_image = convolve(filled_image, kernel)

        # plot
        if save_plot:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

            vmin = 0.1
            vmax2 = np.max(filled_image) * 0.9
            ax1.imshow(image, cmap='Greys', origin='lower', norm=LogNorm(vmin=vmin))
            ax2.imshow(filled_image, cmap='Greys', origin='lower', norm=LogNorm(vmin=vmin, vmax=vmax2))

            # display the plot
            plt.tight_layout()

            # save as a pdf
            plot_name = os.path.join(str(self.name), str(self.name) + '_' +
                                     str(self.wavelength) + '_source_removal.pdf')
            try:
                os.remove(plot_name)
            except:
                pass

            fig.savefig(plot_name, dpi=200)

        return filled_image

    def cuts(self, wavelength, sensor, MCSNR):
        """
        cuts only downloads wanted images from IRAC and removes unwanted images from MIPS
        
        Input
        ==========================
        MCSNR (str):
            name of SNR
            
        sensor (str):
            sensor name
            
        wavelength (str):
            waveband for sensor
            
        """

        sensor_l=len(sensor)+1
        sens=sensor.rjust(sensor_l)
        name = str(sens) + ' ' + str(wavelength) + 'um'
        my_image_files = glob.glob(os.path.join(str(self.name) + '_' + str(wavelength), '*.fits'))
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
            if not os.path.exists(os.path.join(MCSNR)):
                os.makedirs(os.path.join(MCSNR))
            if os.path.exists(os.path.join(my_image_files[i])):
                copyfile(os.path.join(my_image_files[i]), os.path.join(MCSNR + '/' + filename))

        rmtree(os.path.join(str(self.name) + '_' + str(wavelength) + '/'))

    def plot_image(self, data, wcs, coords=None):
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
        vmin = 1e-1
        vmax = np.max(data) * 0.9

        # plot
        axs.imshow(data, cmap='magma', interpolation='nearest', origin='lower', norm=LogNorm(vmin=vmin, vmax=vmax))

        #axs.scatter(coords[0], coords[1], transform=axs.get_transform('fk5'), s=20000, lw=2,
        #            edgecolor='white', facecolor='none')

        # define the apertures to plot (these are the same as in other functions)
        # SNR
        position = SkyCoord(self.ra_deg * u.degree, self.dec_deg * u.degree, frame='icrs')
        snr_aperture = SkyCircularAperture(position, r=self.radius * u.arcsec)
        snr_pix_aperture = snr_aperture.to_pixel(wcs)
        snr_pix_aperture.plot(color='white', lw=5, alpha=1.0)

        # BKG
        r = self.radius * u.arcsec
        r_in = r + (40. * u.arcsec)
        r_out = r + (100. * u.arcsec)
        bkg_ap = SkyCircularAnnulus(position, r_in=r_in, r_out=r_out)
        bkg_ap_pix = bkg_ap.to_pixel(wcs)
        bkg_ap_pix.plot(color='red', lw=3, ls='--', alpha=0.9)


        # MIRI Imager
        r = self.radius * u.arcsec
        r_in = r  - (37. * u.arcsec)
        r_out = r  + (37. * u.arcsec)
        arc_ap = SkyCircularAnnulus(position, r_in=r_in, r_out=r_out)
        arc_ap_pix = arc_ap.to_pixel(wcs)
        arc_ap_pix.plot(color='cyan', lw=3, ls=':', alpha=0.9)

        axs.set_facecolor('black')
        axs.coords.grid(True, color='white', ls='dotted')
        axs.coords[0].set_axislabel('Right Ascension (J2000)')
        axs.coords[1].set_axislabel('Declination (J2000)')

        # indicative MIRI FOV
        x = Angle(self.radius, u.arcsec)
        y = x.degree
        t3_pos = SkyCoord((self.ra_deg) * u.degree, (self.dec_deg + y) * u.degree, frame='icrs')
        ap3 = SkyRectangularAperture(t3_pos, 74 * u.arcsec, 113 * u.arcsec, 0 * u.degree)
        ap3_pix = ap3.to_pixel(wcs)
        ap3_pix.plot(color='yellow', lw=3, ls='-', alpha=0.9)

        axs.set_facecolor('black')
        axs.coords.grid(True, color='white', ls='dotted')
        axs.coords[0].set_axislabel('Right Ascension (J2000)')
        axs.coords[1].set_axislabel('Declination (J2000)')

        # display the plot
        plt.tight_layout()

        # save as a pdf
        plot_name = os.path.join(str(self.name), str(self.name) + '_' +
                                 str(self.wavelength) + '.pdf')
        try:
            os.remove(plot_name)
        except:
            pass

        fig.savefig(plot_name, dpi=200)

    def show_image(self, MCSNR, test_src_coord, wavelength):
        """
        show_image plots cut images for the user to see when running the program
        """
        my_image_file = os.path.join(str(self.name), str(self.sensor) + '_' + str(self.wavelength) + '.fits')

        # load the file data, header, and wcs
        with fits.open(my_image_file) as hdulist:
            header = hdulist[0].header
            data = hdulist[0].data
            wcs = WCS(header)
            size = self.radius * 4
            arc_pix = header['PXSCAL1']
            if arc_pix < 0:
                arc_pix = -arc_pix
            size_pix = size / arc_pix
            factor = 1
            if size_pix in range(100,200):
                factor = (2160/162)
            if size_pix in range(300,400):
                factor = (2160/324)
            if size_pix in range(500,600):
                factor = (2160/529)

            # changing to two twice the SNR diameter
            im_size = (self.radius * 4) / arc_pix
            im_size = int(im_size)
               
            position = SkyCoord(test_src_coord[0] * u.degree, test_src_coord[1] * u.degree, frame='icrs')
            cutout = Cutout2D(data, position, (im_size), wcs=wcs)

            # save the cutout to fits
            fits_name = os.path.join(str(self.name), str(self.name) + '_' +
                                     str(self.wavelength) + '_cut.fits')
            new_hdu = fits.PrimaryHDU(data=cutout.data)
            new_hdulist = fits.HDUList([new_hdu])
            new_hdulist.writeto(fits_name, overwrite=True)

            # fill the point sources
            cutout.data = self.fill_sources(cutout.data, save_plot=True)

            # save the filled image to fits
            fits_name = os.path.join(str(self.name), str(self.name) + '_' +
                                     str(self.wavelength) + '_cut_ps-removed.fits')
            new_hdu = fits.PrimaryHDU(data=cutout.data)
            new_hdulist = fits.HDUList([new_hdu])
            new_hdulist.writeto(fits_name, overwrite=True)

            self.plot_image(cutout.data, cutout.wcs, coords=test_src_coord)

        return my_image_file

    def ap_phot(self, my_image_file, Rad, test_src_coord, wavelength):
        """
        ap_phot applies aperture photometry to calculate the flux of the source. It creates a circle aperture around
        the source referencing the radius of the snr as the radius of the circular aperture. It also creates 4 more
        apertures away from the source which are then averaged and subtracted from the target flux to remove background

        Input
        ==========================
        my_image_files (.fits):
            .fits file from coordinates

        Rad (float):
            radius in units of arcsec

        test_src_coord (float):
            coordinates of target in decimal degrees

        wavelength (float):
            waveband in units of um

        Output
        ==========================
        FluxIR (float):
            calculated flux in units of erg/s/cm^2
        """

        fluxes = []

        r = Rad * u.arcsec
        r_in = r + (40. * u.arcsec)
        r_out = r + (100. * u.arcsec)

        # load the file data, header, and wcs
        with fits.open(my_image_file) as hdulist:
            my_hdu = hdulist[0]
            my_hdu.data = np.nan_to_num(my_hdu.data)
            pixel_area = my_hdu.header['PXSCAL1']**2

            print('Running aperture photometry {}-{} um'.format(self.name, str(wavelength)))
            position = SkyCoord(test_src_coord[0] * u.degree, test_src_coord[1] * u.degree, frame='icrs')
            apertures = SkyCircularAperture(position, r=Rad * u.arcsec)
            phot_table = aperture_photometry(my_hdu, apertures)
            fluxes.append(phot_table['aperture_sum'])
            #print(phot_table['aperture_sum'])


            print('Running background aperture photometry {}-{} um'.format(self.name, str(wavelength)))

            bkg_ap = SkyCircularAnnulus(position, r_in=r_in, r_out=r_out)
            phot_table_bkg = aperture_photometry(my_hdu, bkg_ap)
            #print(phot_table_bkg['aperture_sum'])

            # bkg subtract
            flux_bkg_sub = phot_table['aperture_sum'] - phot_table_bkg['aperture_sum']

        # convert to Jy
        if flux_bkg_sub.value > 0:
            ujy_arcsec = flux_bkg_sub.value * 23.5045
            Jy = ujy_arcsec * pixel_area * 1.e-06
            Jy = Jy[0]

            erg = (Jy * 1.e-23) * (2.997924e14 / ((wavelength)**2))

            self.flux_Jy = Jy
            self.flux_erg = erg

        else:
            erg = 0.0

            self.flux_Jy = 0.0
            self.flux_erg = 0.0

        return erg
    
    def flux_arc(self, my_image_file, Rad, test_src_coord, wavelength):

        r = Rad * u.arcsec
        arc_r_in = r - (37. * u.arcsec)
        arc_r_out = r + (37. * u.arcsec)
        bkg_r_in = r * 1.1
        bkg_r_out = r * 1.25

        # load the file data, header, and wcs
        with fits.open(my_image_file) as hdulist:
            my_hdu = hdulist[0]
            my_hdu.data = np.nan_to_num(my_hdu.data)
            pixel_area = my_hdu.header['PXSCAL1'] ** 2

            position = SkyCoord(test_src_coord[0] * u.degree, test_src_coord[1] * u.degree, frame='icrs')

            print('Running arc aperture photometry {}-{} um'.format(self.name, str(wavelength)))
            arc_ap = SkyCircularAnnulus(position, r_in=arc_r_in, r_out=arc_r_out)
            phot_table_arc = aperture_photometry(my_hdu, arc_ap)

            bkg_ap = SkyCircularAnnulus(position, r_in=bkg_r_in, r_out=bkg_r_out)
            phot_table_bkg = aperture_photometry(my_hdu, bkg_ap)

            # bkg subtract
            flux_bkg_sub = phot_table_arc['aperture_sum'] - phot_table_bkg['aperture_sum']

        # convert to Jy
        if flux_bkg_sub.value > 0:

            ujy_arcsec = flux_bkg_sub.value * 23.5045
            Jy = ujy_arcsec * pixel_area * 1.e-06
            Jy = Jy[0]

            erg = (Jy * 1.e-23) * (2.997924e14 / ((wavelength)**2))

            # we must divide the circumference by the MIRI imager FOV height to get
            # average flux in a MIRI Imager FOV
            circum = 2 * np.pi * Rad
            n_miri_ima = circum / 113.

            self.arc_flux_Jy = Jy / n_miri_ima
            self.arc_flux_erg = erg / n_miri_ima

        else:
            self.arc_flux_Jy = 0.0
            self.arc_flux_erg = 0.0

        return self.arc_flux_erg
        

    def run(self):
        """
        run the class methods
        """
        # Eq2Deg Function
        test_src_coord = self.eq2deg(self.ra, self.dec)

        # Query Function
        self.query(test_src_coord, self.sensor, self.wavelength)

        # Cuts Function
        self.cuts(self.wavelength, self.sensor, self.name)

        # Plot Image function
        my_image_file = self.show_image(self.name, test_src_coord, self.wavelength)

        # Ap Phot Function
        flux = self.ap_phot(my_image_file, self.radius, test_src_coord, self.wavelength)
        
        flux_arc = self.flux_arc(my_image_file, self.radius, test_src_coord, self.wavelength)
        
        return self.flux_Jy, self.flux_erg, self.arc_flux_Jy, self.arc_flux_erg
