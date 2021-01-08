#remove parameterize keyword from all lines in FF
#
#

from openforcefield.typing.engines.smirnoff import ForceField
ff = ForceField('param_valence.offxml', allow_cosmetic_attributes=True)
for handler in ff.registered_parameter_handlers:
    print(handler)
    if len(ff[handler].parameters) == 0:
        continue
    print(ff[handler].parameters)
    for par in ff[handler].parameters:
        if hasattr(par, '_parameterize'):
            par.delete_cosmetic_attribute('parameterize')
        #if hasattr(ff[handler].parameters[0], '_parameterize'):
        #ff[handler].parameters[0].delete_cosmetic_attribute('parameterize')
print(ff)

#<Proper smirks="[*:1]~[#6X3:2]-[#7X3:3]-[#6X4,#1:4]" periodicity1="2" phase1="180.0 * degree" k1_bondorder1="1*kilocalories_per_mole" k1_bondorder2="10*kilocalories_per_mole" id="TIG9" idivf1="1.0"></Proper>

torsionhandler = ff['ProperTorsions']
torsionhandler.add_parameter({'smirks':'[*:1]~[*:2]~[*:3]~[*:4]', 'periodicity1':2, 'phase1':'180.0 * degree', 'k1_bondorder1':'1 * kilocalories_per_mole',  'k1_bondorder2':'1 * kilocalories_per_mole', 'id':'TIG0', 'idivf1':1})
torsionhandler.parameters[-1].add_cosmetic_attribute('parameterize', 'k1_bondorder1,k1_bondorder2')

ff.to_file('test.offxml')
