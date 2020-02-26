import numpy as np
from PIL import Image
import sys, os
sys.path.append('/home/vlb9398/rsz/utilities')
import config
import matplotlib.pyplot as plt

def save_final_im(sgen,nsim,clusname,testflag=0):

    bands = ['PSW','PMW','PLW']

    for i in range(3): # make one set of figs for each band
        # this is all to just make a stupid white figure so things line up
        test_fig = 'pcat_residuals/'+ clusname + '_resid_' + bands[i] + '_' + nsim + '.png'
        test_im = Image.open(test_fig)
        plt.gca().set_axis_off()
        cmap = plt.cm.OrRd
        cmap.set_under(color='white')
        plt.imshow(np.zeros(test_im.size),cmap=cmap,vmin=0.1,vmax=1.0)
        plt.savefig(config.OUTPUT + 'add_sz/' + clusname + '_ibmodel_' + bands[i] + '_fake.png')
        plt.clf()

        if sgen != None :
            list_im = [
                        'pcat_residuals/'+ clusname + '_resid_' + bands[i] + '_' + nsim + '.png',
                        'pcat_residuals/'+ clusname + '_mask_resid_' + bands[i] + '_' + nsim + '.png',
                        'radial_averages/'+ clusname + '_radav_' + bands[i] + '_' + nsim + '.png',
                        'rings/'+ clusname + '_rings_' + bands[i] + '_' + nsim + '.png',
                        'radial_averages/'+ clusname + '_ra_v_r_' + bands[i] + '_' + nsim + '.png',
                        'IB_fits/'+ clusname + '_dI_fit_' + bands[i] + '_' + nsim + '.png'
                        ]
            if i != 0 :
                list_im.append('sim_sz/'+ clusname + '_sze_' + bands[i] + '_' + nsim + '.png')
                list_im.append('corr_comps/'+ clusname + '_xcomps_' + bands[i] + '_' + nsim + '.png')
                list_im.append('corr_comps/'+ clusname + '_xclean_' + bands[i] + '_' + nsim + '.png')

                if testflag != 0: # there is an extra plot from fake IB model in this case
                    list_im.append('add_sz/' + clusname + '_ibmodel_' + bands[i] + '_fake' + '.png')
                    list_im.append('add_sz/' + clusname + '_ibmodel_' + bands[i] + '_' + nsim + '.png')
                    list_im.append('add_sz/' + clusname + '_ibmodel_' + bands[i] + '_fake' + '.png')

        else :
            list_im = [
                        'pcat_residuals/'+ clusname + '_resid_' + bands[i] + '_real.png',
                        'pcat_residuals/'+ clusname + '_mask_resid_' + bands[i] + '_real.png',
                        'radial_averages/'+ clusname + '_radav_' + bands[i] + '_real.png',
                        'rings/'+ clusname + '_rings_' + bands[i] + '_real.png',
                        'radial_averages/'+ clusname + '_ra_v_r_' + bands[i] + '_real.png',
                        'IB_fits/'+ clusname + '_dI_fit_' + bands[i] + '_real.png'
                        ]
            if i != 0 :
                list_im.append('sim_sz/'+ clusname + '_sze_' + bands[i] + '_real.png')
                list_im.append('corr_comps/'+ clusname + '_xcomps_' + bands[i] + '_real.png')
                list_im.append('corr_comps/'+ clusname + '_xclean_' + bands[i] + '_real.png')

                if testflag != 0: # there is an extra plot from fake IB model in this case
                    list_im.append('add_sz/' + clusname + '_ibmodel_' + bands[i] + '_fake.png')
                    list_im.append('add_sz/' + clusname + '_ibmodel_' + bands[i] + '_real.png')
                    list_im.append('add_sz/' + clusname + '_ibmodel_' + bands[i] + '_fake.png')

        list_im = [config.OUTPUT + s for s in list_im] # add the full file path at the start of each file name

        # create first set of 3 horizontal images
        if len(list_im)%3 != 0: # we only want 3 images at a time
            list_len = round(len(list_im)/3)
        else :
            list_len = len(list_im)/3

        l = 0
        imgs_horz = []
        for j in range(0,len(list_im),3):
            imgs = [Image.open(i) for i in list_im[j:j+3]]

            # pick the image which is the smallest, and resize the others to match it
            min_shape = sorted([(np.sum(i.size), i.size) for i in imgs])[0][1]
            imgs_comb = np.hstack((np.asarray(i.resize(min_shape)) for i in imgs))

            # save that beautiful picture
            imgs_comb = Image.fromarray(imgs_comb)
            imgs_comb.save(config.OUTPUT + 'final/rsz_comp_%s.png'%(l))
            imgs_horz.append(config.OUTPUT + 'final/rsz_comp_%s.png'%(l))
            l += 1

        # vertically stack the 3 horizontal saved images we made
        imgs_h = [Image.open(i) for i in imgs_horz]
        max_shapeh = sorted([(np.sum(i.size), i.size) for i in imgs_h])[-1][1]
        imgs_combh = np.vstack(imgs_h)
        imgs_combh = Image.fromarray(imgs_combh)
        imgs_combh.save(config.OUTPUT + 'final/%s_rsz_final_%s_%s.png' %(clusname,nsim,bands[i]))

        # go through and delete all of the files that we no longer need
        for s in imgs_horz:
            os.remove(s)

        #delete originals
        # for m in list_im :
        #     os.remove(m)

if __name__ == '__main__' :
    sgen = 3
    nsim = 0
    clusname = 'rxj1347'
    save_final_im(str(sgen),str(nsim),clusname,testflag=1)
