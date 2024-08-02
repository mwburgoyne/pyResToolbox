from pyrestoolbox.classes import z_method, c_method, pb_method, rs_method, bo_method, uo_method, deno_method, co_method, kr_family, kr_table, class_dic
import sys

def validate_methods(names, variables):
    for m, method in enumerate(names):
        #print(class_dic['cmethod'][cmethod.upper()])
        if type(variables[m]) == str:
            try:
                variables[m] = class_dic[method][variables[m].upper()]
            except:
                print("An incorrect method was specified")
                sys.exit()
    if len(variables) == 1:
        return variables[0]
    else:
        return variables