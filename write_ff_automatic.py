from mkutils.gromacs_writer import FFWriter
import json, os

ff_file = 'sdk_parameters.json'
in_file = 'spica_par.json'

species = {
           'W': {'weight': 54.045, 'at_num': 30},
           'SOD': {'weight': 54.045, 'at_num': 30},
           'CLA': {'weight': 54.045, 'at_num': 30},
           'SO4': {'weight': 54.045, 'at_num': 30},
           'CM': {'weight': 54.045, 'at_num': 30},
           'CT': {'weight': 54.045, 'at_num': 30}
          }

def K(kcal): 
     R = 8.3144598 
     kcal_to_J = 4186.8 
     return kcal*kcal_to_J/R 
  
if not os.path.isfile(ff_file): 
    with open(ff_file, 'w') as f: 
        _dict = {'atomtypes': {}, 'crossints': {}}                          
        json.dump(_dict, f) 
           
with open(in_file, 'r') as f: 
    spica = json.load(f) 
 
ff = FFWriter(ff_file, rule='SDK') 
 
# self interactions 
all_species = list(species.keys())
for item in all_species: 
    for param in spica.get('params'): 
        if param.get('param') == 'pair': 
            types = param.get('types') 
            if types[0] == item and types[1] == item: 
                name = 'sdk{:s}'.format(item)  
                potential = param.get('potential')[2:].split('_') 
                lr, la = float(potential[0]), float(potential[1]) 
                eps = K(param.get('epsilon')) 
                sig = param.get('sigma')/10. 
                weight = species.get(item).get('weight') 
                at_num = species.get(item).get('at_num') 
                ff.add_atomtype(name, lr, la, eps, sig, weight, at_num) 
 
def condition(types, itemi, itemj):
    if itemj == itemi:
        if itemj == types[0] and itemi == types[1]:
            return True
        else:
            return False
    if itemj in types and itemi in types:
        return True
    else:
        return False
        
for i, itemi in enumerate(all_species):
    for itemj in all_species[i:]:
        for param in spica.get('params'): 
            if param.get('param') == 'pair': 
                types = param.get('types')
                if condition(types, itemi, itemj): 
                    itemi_name = 'sdk{:s}'.format(itemi)  
                    itemj_name = 'sdk{:s}'.format(itemj)  
                    eps = K(param.get('epsilon')) 
                    sig = param.get('sigma')/10.
                    ff.add_crossint(itemi_name, itemj_name, 
                                    eps_mix=eps, sig_mix=sig
                                    )
                    break
        

 
ff.write_tables(shift=True, cutoff=1.5) 
ff.write_forcefield(outfile='ffnonbonded.itp') 
