import numpy as np
import PIL
import sys
sys.path.append('utilities')
import config

def save_final_im(sgen,nsim,clusname,testflag=0):

    bands = ['PSW','PMW','PLW']

    for i in range(3): # make one set of figs for each band
        if sgen != None :
            list_im = [
                        'sim_sz/'+ clusname + '_sze_' + bands[i] + '_' + nsim + '.png',
                        'pcat_residuals/'+ clusname + '_resid_' + bands[i] + '_' + nsim + '.png',
                        'pcat_residuals/'+ clusname + '_mask_resid_' + bands[i] + '_' + nsim + '.png',
                        'corr_comps/'+ clusname + '_xcomps_' + bands[i] + '_' + nsim + '.png',
                        'corr_comps/'+ clusname + '_xclean_' + bands[i] + '_' + nsim + '.png',
                        'radial_averages/'+ clusname + '_radave_' + bands[i] + '_' + nsim + '.png',
                        'rings/'+ clusname + '_rings_' + bands[i] + '_' + nsim + '.png',
                        'radial_averages/'+ clusname + '_ra_v_r_' + bands[i] + '_' + nsim + '.png',
                        'IB_fits/'+ clusname + '_dI_fit_' + bands[i] + '_' + nsim + '.png'
                        ]
            if testflag != 0: # there is an extra plot from fake IB model in this case
                list_im.append('add_sz/' + clusname + '_ibmodel_' + bands[i] + '_' + nsim + '.png')

        else :
            list_im = [
                        'sim_sz/'+ clusname + '_sze_' + bands[i] + '_real.png',
                        'pcat_residuals/'+ clusname + '_resid_' + bands[i] + '_real.png',
                        'pcat_residuals/'+ clusname + '_mask_resid_' + bands[i] + '_real.png',
                        'corr_comps/'+ clusname + '_xcomps_' + bands[i] + '_real.png',
                        'corr_comps/'+ clusname + '_xclean_' + bands[i] + '_real.png',
                        'radial_averages/'+ clusname + '_radave_' + bands[i] + '_real.png',
                        'rings/'+ clusname + '_rings_' + bands[i] + '_real.png',
                        'radial_averages/'+ clusname + '_ra_v_r_' + bands[i] + '_real.png',
                        'IB_fits/'+ clusname + '_dI_fit_' + bands[i] + '_real.png'
                        ]
            if testflag != 0: # there is an extra plot from fake IB model in this case
                list_im.append('add_sz/' + clusname + '_ibmodel_' + bands[i] + '_real.png')

        list_im = [config.OUTPUT + s for s in list_im] # add the full file path at the start of each file name


    '''
        We need a good way of deciding what to do with the last, 10th
        plot in case the testing IB model is used.
        This code will only do the 9 case.
    '''
    # create first set of 3 horizontal images
    for j in range(3):
        imgs = [PIL.Image.open(i) for i in list_im[j,j+3]]
        # pick the image which is the smallest, and resize the others to match it (can be arbitrary image shape here)
        min_shape = sorted([(np.sum(i.size), i.size) for i in imgs])[0][1]
        imgs_comb = np.hstack((np.asarray(i.resize(min_shape)) for i in imgs))

        # save that beautiful picture
        imgs_comb = PIL.Image.fromarray(imgs_comb)
        imgs_comb.save(config.OUTPUT + 'final/rsz_comp_%s.png'%(j))

    # maybe put a check in here to make sure these files get made...
    imgs_horz = [config.OUTPUT + 'final/rsz_comp_0.png',
                    config.OUTPUT + 'final/rsz_comp_1.png',
                    config.OUTPUT + 'final/rsz_comp_2.png'
                ]
    # vertically stack the 3 horizontal saved images we made
    imgs_h = [PIL.Image.open(i) for i in imgs_horz]
    min_shapeh = sorted([(np.sum(i.size), i.size) for i in imgs_h])[0][1]
    imgs_combh = np.vstack((np.asarray(i.resize(min_shape)) for i in imgs_h))
    imgs_combh = PIL.Image.fromarray(imgs_combh)
    imgs_combh.save(config.OUTPUT + 'final/rsz_comp_final.pdf')

    # go through and delete all of the files that we no longer need
    ''' Lets wait on this until we actually test that the final figure comes out okay
    for s in list_im:
        os.remove(list_im[s])
    '''
