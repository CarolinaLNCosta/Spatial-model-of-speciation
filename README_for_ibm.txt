To run the spatially-explicit IBM utilized in Costa et al. 2018:

1) Compile the file 'ibm.f90' in your operating system.
2) Put the executable file and the input.in file in the same directory.
3) Run the executable file 'ibm.exe'.

Input file:

input.in: contain parameters to initialize the simulation:

- T: number of generations
- M: size of initial population 
- L1 and L2: size of lattice 
- mut: mutation rate
- diff: probability of dispersion 
- Q: probability of individual do not reproduce
- deltat: time interval in which species are calculated
- S: radius defining the mating range
- G: genetic threshold
- B: genome size
- P: number of potential mates inside mating range
- IREAD: running initial or continuing previous simulation 
- njump: number of neighbouring sites used for dispersion, ranging from 4 to 32. 
- nrmax: maximum increment in the radius size if an individual do not find potential mates in this mating range. 

Output files:

- pop.dat: spatial coordinates and genome of all individuals.
- species-time.dat: number of species in function of time.
- species-plot: spatial coordinates and species that each individual belongs. 
- abund.dat: number of individuals in each species at the end of simulation.
- mrcat-phy.dat: MRCAT matrix selecting one individual per species.
- phy-genetics: genome of individuals in mrcat-phy.dat file. 
- mrcat.dat: MRCAT matrix for all individuals at the end of simulation.
