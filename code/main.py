# Loads csv file
#main
line = np.genfromtxt(open('/home/murphyj/Desktop/Coding/SNR_list.csv', "r"), names=True, delimiter=',', dtype=None)

SNR = 34 # position in list
instrument_list = ['IRAC', 'MIPS']
band_list = [3.6, 4.5, 5.8, 8.0, 24, 70, 160]

sensor = instrument_list[0]
wavelength = band_list[1]

a = Angle(im_ra[SNR], u.hour)
im_ra_deg = a.degree

b = Angle(im_dec[SNR])
im_dec_deg = b.degree

test_src_coord = [im_ra_deg, im_dec_deg]

my_image_files = glob.glob(os.path.join(MCSNR[SNR].decode("utf-8"), '*.fits'))
n_images = len(my_image_files)

# cycle through images
for i in my_image_files:
    # load the file data, header, and wcs
    with fits.open(i) as hdulist:
        header = hdulist[0].header
        data = hdulist[0].data
        wcs = WCS(header)
        plot_image(data, wcs, coords=test_src_coord)

for i in my_image_files:
    # load the file data, header, and wcs
    with fits.open(i) as hdulist:
        size = Rad[SNR] * 2
        arc_pix = header['PXSCAL1']
        if arc_pix < 0:
            arc_pix = -arc_pix
        size_pix = size / arc_pix
        header = hdulist[0].header
        data = hdulist[0].data
        wcs = WCS(header)
        position = SkyCoord(test_src_coord[0] * u.degree, test_src_coord[1] * u.degree, frame='icrs')
        cutout = Cutout2D(data, position, (size), wcs=wcs)

        # plt.imshow(cutout.data, origin='lower')
        # plt.tight_layout()
        # plt.show()
        plot_image(cutout.data, cutout.wcs, coords=test_src_coord)

