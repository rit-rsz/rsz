import sys
sys.path.append('source_handling/')
#this should be fine but you may need to change the path
#to whereever get data is in your file system
from clus_get_data import clus_read_file


dir = '/data/mercado/SPIRE/bethermin_sims/a0370/' #change this to whatever your directory for keeping files is
filename = ['a0370_PSW_sim0110.fits', 'a0370_PMW_sim0110.fits', 'a0370_PLW_sim0110.fits'] #this is a list of the filenames.
#please put the file names in order of PSW, PMW, PLW
band = ['PSW', 'PMW', 'PLW']
maps = [[],[],[]]
for i in range(3):
    maps[i] = clus_read_file(dir + filename[i], band[i], 'some cluster name')
#band name should be 'PSW, PMW, PLW'
#clusname doesn't matter unless your script cares about clusname in which cause put
#the name of the cluster you are using,
