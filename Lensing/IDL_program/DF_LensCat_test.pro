@DF_LensCat_FP.pro

FUNCTION testlensing, t

input_cat = 'test_dat_files/test_catalogue_PSW_sim15.dat'
DF_X = '/home/vaughan/rsz/Lensing/IDL_program/alpha_x_ALL_rxj1347_z1p75.txt'
DF_Y = '/home/vaughan/rsz/Lensing/IDL_program/alpha_y_ALL_rxj1347_z1p75.txt'
lens_z = 0.451 ;from csv file jack sent
LL = 150
pixscale = 18
searchedlensed = 0.75
npixels = 290
x = DF_LensCat_FP(input_cat, DF_X, DF_Y, lens_z, pixscale, LL, NPIXELS=npixels)
print, x
END
