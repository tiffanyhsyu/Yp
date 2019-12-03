# yMCMC

yMCMC is a code that solves for the best fit parameters that describe an emission line spectrum, given their
 hydrogen and helium flux ratios. There are 8 parameters yMCMC solves for:
- y+, singly ionized helium abundance
- T_e, the electron temperature [K]
- log10(n_e), the electron density [cm^-3]
- c(Hb), the reddening parameter
- a_H, the amount of underlying hydrogen stellar absorption [A], relative to the amount at Hb
- a_He, the amount of underlying helium stellar absorption [A], relative to the amount at HeI 4026
- tau_He, the helium optical depth parameter, normalized to the value at HeI 3889
- log10(xi), ratio of neutral to singly ionized hydrogen densities


yMCMC steps through this parameter space, predicts the would-be emission line flux ratios given the parameters, and
calculates the log-likelihood function of the model.

The following are instructions on how to set up and run yMCMC on your emission line systems. These steps require
the user to have existing measurements of the emission line flux ratios and equivalent widths (EW) of as many of the
following hydrogen and helium emission line ratios as possible:

- Hydrogen: H8, H-delta, H-gamma, H-beta, H-alpha, H-alpha, P-gamma
- Helium: HeI3889, HeI4026, HeI4472, HeI5016, HeI5876, HeI6678, HeI7065, HeI10830

1. Create an input file per system that lists the wavelength, flux ratio, flux ratio error, EW, EW error for the
measured emission line ratios. The flux ratios need to be normalized to the value at H-beta for optical lines,
and normalized to the value at P-gamma for the NIR line HeI10830. See the files in folder ~/yMCMC/test_data/test/
for examples. Note that there is a +/-3 Angstrom tolerance on the wavelengths listed in the sample files, and the
'Species' column is not required. Create a folder that these input files are placed in, for example,
'test/data/galsamp/'.

2. In the galaxy.py file, create a new function that will read in the input file(s) you have from Step 1. For a
galaxy sample named 'galsamp' with two systems named 'galaxy1' and 'galaxy2', the function would look like:
```
    def load_galsamp(galaxyname):
        outdict = dict()
        dir = '/test_data/galsamp/'
        if galaxyname == 'galaxy1':
            gname = 'gal1_input'
            T_OIII = 10000.
            full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
        elif galaxyname == 'galaxy2':
            gname = 'gal2_input'
            T_OIII = 10000.
            full_tbl = Table.read(os.getcwd() + dir + gname, format='ascii', delimiter=' ')
        else:
            print('Galaxy not known: {0:s}'.format(galaxyname))
            return None
        outdict['full_tbl'] = full_tbl
        outdict['T_OIII'] = T_OIII
        return outdict
```
The number of if statements must match the number of systems in your sample. Here, 'gname' refers to the name of the
input file created in Step 1 for your galaxy. Make sure that these input files are placed in the '/test_data/galsamp'
folder, or that you edit the path in the directory line to match where the files are stored. 'T_OIII' is the
temperature prior for that galaxy, which the user should measure via the direct method and using the [OIII] lines.
The values 10000. given here are placeholders.

3. Go to mcmc_all.py and edit the following line (near the top of the file, Line ~14) such that the function
created in Step 2 is called on:
```
    galdict = galaxy.load_galsamp(self.galaxyname)
```
4. Also in mcmc_all.py, towards the middle (Line ~320), add the names of your galaxies in galsamp to a new list:
```
    galsamp = ['galaxy1', 'galaxy2']
```
Directly under the lists, add/change the following line to read in your new list:
```
    rungal = 'galsamp'
```
5. Also in mcmc_all.py, ~10 lines below Step 4, add your galaxy sample to the series of if statements:
```
    elif rungal == 'galsamp':
        # Run MCMC on our galaxy sample
        galfail = []
        for gal in galsamp:
            try:
                print ('Working on', gal)
                MCMCgal(gal)
            except IOError:
                print('ERROR :: The following galaxy data could not be found: {0:s}'.format(gal))
                galfail += [gal]
            except TypeError:
                print('ERROR :: The following galaxy is not known: {0:s}'.format(gal))
                galfail += [gal]
            except ValueError:
                print('ERROR :: The following galaxy failed: {0:s}'.format(gal))
                galfail += [gal]
        print('The following galaxies failed:\n' + '\n'.join(galfail))
```
6. Once these additions/changes are made, you can run mcmc_all.py. The MCMC by default uses 500 walkers, each
taking 1000 steps. The MCMC chains will be saved in files ending with '_500walkers_1000steps.npy'.