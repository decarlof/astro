from astropy.io import fits


fits_image_filename = '/Users/decarlo/conda/astro/data//hmi.meanpf_720s.20170108_192400_TAI.mf_br.fits'
hdul = fits.open(fits_image_filename)

print('**********************')
print(hdul)
print('**********************')
print(hdul.info())
print('**********************')

print(hdul[0].data)

a = hdul[0].data

print(a*2)

# fit_header = hdul[0].header
# for i in range(0, len(fit_header), 1):
#     header = hdul[0].header[i]
#     print(i, header)

# fit_header = hdul[0].header
# for key in fit_header:
#     print(key, fit_header[key])


# print(fit_header['CHECKSUM'])
# print(fit_header[7])


hdul.close()


# preferred
with fits.open(fits_image_filename) as hdul:
    print(hdul[0].data)
