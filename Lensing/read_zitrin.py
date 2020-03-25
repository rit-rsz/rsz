import numpy as np

def read_zitrin(npix,file):

    z_array = []

    with open(file,'r') as f :
        data = f.readlines()
        
        for i in range(len(data)) :
            z_array.append([float(x) for x in data[i].strip(' ').split()][:-1])

    f.close()

    z_array = np.resize(np.array(z_array),(npix[0],npix[0]))

    return z_array
