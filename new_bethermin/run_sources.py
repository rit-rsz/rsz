from bethermin12_sim import gencat
gn = gencat()
ngen = 10000 # Generate 10,000 sources
cat = gn.generate(ngen)  # Generates sources, but not flux densities
wavearr = [250.0, 350.0, 500.0] #In um
cat = gn.generate(ngen, wave=wavearr) # Generates sources and flux densities
