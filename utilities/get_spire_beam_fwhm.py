

def get_spire_beam_fwhm(band):
    if band == 'PSW':
        beamFWHM = 18.0
        return beamFWHM
    elif band == 'PMW':
        beamFWHM = 25.0
        return beamFWHM
    elif band == 'PLW':
        beamFWHM = 36.0
        return beamFWHM
    elif band == 'BOLOCAM':
        beamFWHM = 144.0
        return beamFWHM
    else:
        print('Unknown band:', band)
