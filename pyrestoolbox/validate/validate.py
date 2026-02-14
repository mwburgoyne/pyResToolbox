from pyrestoolbox.classes import z_method, c_method, pb_method, rs_method, bo_method, uo_method, deno_method, co_method, kr_family, kr_table, class_dic

        
def validate_methods(names, variables):
    # Backward compatibility mapping
    legacy_map = {'BUR': 'BNS'}
    
    for m, method in enumerate(names):
        #print(class_dic['cmethod'][cmethod.upper()])
        if isinstance(variables[m], str):
            # Replace legacy method names with current ones
            method_str = variables[m].upper()
            method_str = legacy_map.get(method_str, method_str)

            try:
                variables[m] = class_dic[method][method_str]
            except KeyError:
                raise ValueError(f"An incorrect method was specified: '{variables[m]}' for '{method}'")
    if len(variables) == 1:
        return variables[0]
    else:
        return variables