# utilities
#main
import os, glob
from pprint import pprint
from shutil import copyfile, rmtree
from photutils import SkyCircularAperture
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

        return im_ra_deg, im_dec_deg

    def query(self, test_src_coord, sensor, wavelength):
        """
        query creates a string to reference the image database
        
        Input
        ==========================
        test_src_coord (float):
            coordinates of SNR in decimal degrees
            
        sensor (str):
            sensor name
            
        wavelength (str):
            waveband for sensor
            
        Output
        ==========================
        querys SHA database and downloads .fits files into test folders
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
        my_image_files = glob.glob(os.path.join(str(self.name) + '_' + str(wavelength), '*.fits'))
        counter = 0
        filename = sens[1:] + '_' + str(wavelength) + '.fits'
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
        my_image_files = glob.glob(os.path.join(MCSNR, '*.fits'))
        for f in my_image_files:
            if str(wavelength) in f:
                my_image_file = f

        # load the file data, header, and wcs
        with fits.open(my_image_file) as hdulist:
            header = hdulist[0].header
            data = hdulist[0].data
            wcs = WCS(header)
            size = self.radius * 4
               
            position = SkyCoord(test_src_coord[0] * u.degree, test_src_coord[1] * u.degree, frame='icrs')
            cutout = Cutout2D(data, position, (size), wcs=wcs)

            self.plot_image(cutout.data, cutout.wcs, coords=test_src_coord)
        return my_image_files

    def ap_phot(self, my_image_files, Rad, test_src_coord, wavelength):
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
        fluxesav = []
        for i in my_image_files:
    
            # load the file data, header, and wcs
            with fits.open(i) as hdulist:
                my_hdu = hdulist[0]
                my_hdu.data = np.nan_to_num(my_hdu.data)
        
                position = SkyCoord(test_src_coord[0] * u.degree, test_src_coord[1] * u.degree, frame='icrs')
                apertures = SkyCircularAperture(position, r=Rad * u.arcsec)

                phot_table = aperture_photometry(my_hdu, apertures)

                fluxes.append(phot_table['aperture_sum'])
        
                x = Angle(Rad, u.arcsec)
                y = 2 * x.degree
        
                t1_pos = SkyCoord((test_src_coord[0] + y) * u.degree, (test_src_coord[1] + y) * u.degree, frame='icrs')
                t2_pos = SkyCoord((test_src_coord[0] + y) * u.degree, (test_src_coord[1] - y) * u.degree, frame='icrs')
                t3_pos = SkyCoord((test_src_coord[0] - y) * u.degree, (test_src_coord[1] + y) * u.degree, frame='icrs')
                t4_pos = SkyCoord((test_src_coord[0] - y) * u.degree, (test_src_coord[1] - y) * u.degree, frame='icrs')
        
                ap1 = SkyCircularAperture(t1_pos, r=Rad * u.arcsec)
                ap2 = SkyCircularAperture(t2_pos, r=Rad * u.arcsec)
                ap3 = SkyCircularAperture(t3_pos, r=Rad * u.arcsec)
                ap4 = SkyCircularAperture(t4_pos, r=Rad * u.arcsec)
        
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

        unit = flux * 23.5045
        ap_ar = np.pi * Rad**2
        jan = unit / ap_ar
        erg = jan * 10**-17

        if wavelength == 3.6:
            band = 9.671e13 - 7.687e13
        if wavelength == 4.5:
            band = 7.687e13 - 5.996e13
        if wavelength == 5.8:
            band = 6.118e13 - 4.835e13
        if wavelength == 8.0:
            band = 4.835e13 - 3.224e13
        if wavelength == 24:
            band = 1.499e13 - 1.070e13
        if wavelength == 70:
            band = 5.996e12 - 3.331e12
        if wavelength == 160:
            band = 2.306e12 - 1.578e12
        fluxIR = erg * band
        FluxIR = fluxIR[0]

        return FluxIR

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
        my_image_file = self.show_image(self.radius, self.name, test_src_coord, self.wavelength)

        # Ap Phot Function
        flux = self.ap_phot(my_image_file, self.radius, test_src_coord, self.wavelength)
        
        return flux[0]
