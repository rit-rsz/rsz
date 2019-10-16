bethermin12_sim
=============================

A package to generate realizations of the Bethermin et al. 2012
far-IR phenomenological model.

###Installation
The usual

	python setup.py install

###Usage
To generate sources:

	from bethermin12_sim import gencat
	gn = gencat()
	ngen = 10000 # Generate 10,000 sources
	cat = gn.generate(ngen)  # Generates sources, but not flux densities
	wavearr = [250.0, 350.0, 500.0] #In um
	cat = gn.generate(ngen, wave=wavearr) # Generates sources and flux densities

To generate simulated maps using a Gaussian beam:

	from bethermin12_sim import genmap_gauss
	gm = genmap_gauss()
	area = 0.25
	maps = gm.generate(area, verbose=True)

Further information is available through the usual python help mechanism.

### Dependencies
This depends on a number of python packages:
* [numpy](http://www.numpy.org/)
* [scipy](http://www.scipy.org/)
* [astropy](http://www.astropy.org/)

### References
* The Bethermin 2012 model is described in
  [Bethermin et al. (2012)](http://dx.doi.org/10.1088/2041-8205/757/2/L23)




