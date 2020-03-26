@DF_LensCat_FP.pro

FUNCTION testboy, t

input_cat = 'Input_params/sides_PLW_sim.npy'
DF_X = 'Input_Params/alpha_x_ALL_rxj1347_z1p75.txt'
DF_Y = 'Input_Params/alpha_y_ALL_rxj1347_z1p75.txt'
lens_z = 0.451 ;from csv file jack sent
LL = 150
pixscale = 18
searchedlensed = 0.75
npixels = 22500

x = DF_LensCat_FP(input_cat, DF_X, DF_Y, lens_z, pixscale, LL)
print, x
END
