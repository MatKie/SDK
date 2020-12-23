from mkutils.gromacs_writer import FFWriter
import json, os

ff_file = 'sdk_parameters.json'

def K(kcal):
    R = 8.3144598
    kcal_to_J = 4186.8
    return kcal*kcal_to_J/R

if not os.path.isfile(ff_file):
    with open(ff_file, 'w') as f:
        _dict = {'atomtypes': {}, 'crossints': {}}       
        json.dump(_dict, f)

ff = FFWriter(ff_file, rule='SDK')

ff.add_atomtype('sdkW', 12., 4., K(0.7), 0.4371, 54.045, 30)
ff.add_atomtype('sdkSOD', 12., 4., K(0.35), 0.4371, 77.035, 41)
ff.add_atomtype('sdkCLA', 12., 4., K(0.35), 0.4371, 71.483, 37)
ff.add_atomtype('sdkSO4', 9., 6., K(0.7), 0.4321, 96.06 , 48)
ff.add_atomtype('sdkCM', 9., 6., K(0.42), 0.4506, 42.07914 , 24)
ff.add_atomtype('sdkCT', 9., 6., K(0.469), 0.4585, 43.07914 , 25)


ff.add_crossint('sdkW', 'sdkSOD', eps_mix=K(0.895))
ff.add_crossint('sdkW', 'sdkCLA', eps_mix=K(0.895))
ff.add_crossint('sdkW', 'sdkSO4', eps_mix=K(1.1))
ff.add_crossint('sdkW', 'sdkCM', eps_mix=K(0.34))
ff.add_crossint('sdkW', 'sdkCT', eps_mix=K(0.36))

ff.add_crossint('sdkSOD', 'sdkCLA', eps_mix=K(0.895))
ff.add_crossint('sdkSOD', 'sdkSO4', eps_mix=K(1.1))
ff.add_crossint('sdkSOD', 'sdkCM', eps_mix=K(0.895))
ff.add_crossint('sdkSOD', 'sdkCT', eps_mix=K(0.895))

ff.write_tables(shift=True, cutoff=1.5)
ff.write_forcefield(outfile='ffnonbonded.itp')


