
################################################################################
# NAME : clus_sim_background.py
# DATE STARTED : September 27, 2019
# AUTHORS : Dale Mercado
# PURPOSE : The wrapper for Conley's Bethermin simmulation generator. fluxcut is in Jy(/beam if you like).
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
#
################################################################################

def clus_sim_background(clusname,simradave,genbethermin=1,genpowerlaw=0,fluxcut=0,saveplots=1,savemaps=0,
            genradave=1,addnoise=0,yeslens=1,resolution='fr',nsim=0,bolocam=0,verbose=1,errmsg=None):

        # Welcome the user
        print(' ')
        if verbose:
            'Welcome to background simulator v 1.0 Python Edition'
        print(' ')

        # Set up ps plotting if appropriate
        if saveplots:
            '''Need to create Nuplot for this aswell'''

        # Set up some path handling
        bethfile = str('bethermin_'+ clusname)
        bethcat = str('bethermin_'+ clusname + '_cat.fits')

        # First get the cluster parameters
        if verbose:
            print('Fetching cluster parameters')
        # Use the version of clusparams that we have already worked with
        params, err = clus_get_clusparams(clusname,verbose=verbose)
        if err:
            if verbose:
                print('clus_get_clusparams exited with error: ' + err)
            exit()

        # now get the spire maps
        if verbose:
            print('Fetching SPIRE maps.')
        maps, err = clus_get_data(clusname=clusname, resolution=resolution, bolocam=bolocam, verbose=verbose)
        if err:
            if verbose:
                print('clus_get_data exited with error: ' + err)
            exit()

        #  Now generate the beam maps for the three colors
        if verbose:
            print('Generating beam maps.')

        ncols =  len(maps)
        wave = [250.,350.,500.] # In units of um

        if resolution == 'fr':
            pixsize = np.repeat(2.,ncols)

        else:
            pixsize = [6.,6.,6.]

        if bolocam:
            pixsize = [pixsize,6]
            wave = [wave,2100.]

        # these were string arrays
        psffiles =  []
        bands = []

        for icol in range(ncol):
            if icol < 3:
                bands[icol] = maps['band']
            else:
                bands[icol] = 'BOLOCAM'
            psffiles[icol] = str(CLUSBOX+bands[icol]+'_psf.fits')

            '''This is a new script to look at'''
            SMAP_WRITE_FITSBEAM

        print(' ')

        if genbethermin:
            # Now fetch a bunch of lookup stuff and generate the source
            #  matrix if requested
            if verbose:
                print('Setting up bethermin model.')

            # This fetches the cold andstarburst templates
            tpldir = CLUSDATA + 'lookup/'
            '''another one'''
            cold = MRDFITS(tpldir+'ias_cold.fits',1)
            starburst = MRDFITS(tpldir+'ias_starburst.fits',1)

            # this is required for later if required
            # Need to find out which are inputs or outputs and order of the returns
            if 0: # if what is 0???
                READCOL,'pk500mu_nfw_3Rvir_v230709.out',k_theta_m,cl_m_lensed,$
                cl_m_1h, cl_m_2h, cl_m_tot

            if gendlots:
                # There was code here but it is commented out in the idl version
                # May be depreciated
                something = None
            else:
                dlots = MRDFITS(str(tpldir,'dNdLogLdzdOmega.fits'),1)

            print('Staring Bethermin.')

            # These inputs will have to be changed inorder to account for the pointer things
            # Test is something like the maps for racen
            SMAP_GENERATE_BETHERMIN(dlots,cold,starburst,psffiles,bethfile,AREA=0.25,
                BANDS=bands,WAVE=wave,RACEN=maps[0]['astr']['crval'][0],DECCEN=(*maps[0]).astr.crval[1],
                CATFILE=STRING(!CLUSSBOX,bethcat),OUTDIR=!CLUSSBOX,VERBOSE=verbose)

            if genpowerlaw:
                if verbose:
                    print('Setting up power law model.')

                catfile = 'powerlawcat.fits'
                outfile = 'powerlawimage'

                knotpos = [1e-1,2,5,10,20,45,100,200,450,1000] * 1e-3
                knotval = [10.07,7.05,6.25,5.919,5.139,4.038,2.596,1.42,0.57,-0.45]

                '''Another new script'''
                # Not sure how often this will even be used
                BROKEN_IMAGE(outfile,knotpos,knotval,AREA=0.01,NFLUXES=1e6,racen=(*maps[0]).astr.crval[0],
                deccen=(*maps[0]).astr.crval[1],CATFILE=catfile,PIXSIZE=1.0)

                catmv = str('powerlaw*'+CLUSBOX)
                # There is a spawn function which I'm not sure how to use
                # Ben or Victoria have more experience
                # SPAWN,STRING('mv ',catmv)

                ncols == 1
                bethcat = catfile
                bethfile = outfile

            for icol in range(ncols):
                '''another one'''
                CLUS_FORMAT_BETHERMIN(icol,STRING(!CLUSSBOX,bethcat),STRING(!CLUSSBOX,bethfile),
                    bands[icol],clusname,pixsize[icol],FLUXCUT=fluxcut,ZZERO=params.z,RETCAT=lozcat)

                if yeslens == 1:

                    print('Starting ', clusname, ' lens run for ', bands[icol], '...')
                
