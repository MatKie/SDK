from mkutils.gromacs_writer import FFWriter
import json, os

ff_file = 'saft_parameters.json'

if not os.path.isfile(ff_file):
    with open(ff_file, 'w') as f:
        _dict = {'atomtypes': {}, 'crossints': {}}       
        json.dump(_dict, f)

ff = FFWriter(ff_file)

ff.add_atomtype('Wift25', 8., 6., 305.21, 0.29016, 18.015028, 10)

ff.add_atomtype('CT', 15.947, 6., 358.37, 0.45012, 43.08698, 25)


ff.add_crossint('Wift25', 'CT', k=0.31)

ff.write_tables(shift=True, cutoff=1.5)
ff.write_forcefield(outfile='ffnonbonded.itp')


