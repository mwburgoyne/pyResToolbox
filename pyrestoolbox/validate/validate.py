"""
Method string-to-Enum validation with legacy name mapping.

Functions
---------
validate_methods    Convert string method specifications to Enum instances
"""

__all__ = ['validate_methods', 'validate_choice']

from pyrestoolbox.classes import z_method, c_method, pb_method, rs_method, bo_method, uo_method, deno_method, co_method, kr_family, kr_table, class_dic


def validate_methods(names, variables):
    # Backward compatibility mapping
    legacy_map = {'BUR': 'BNS'}

    for m, method in enumerate(names):
        if isinstance(variables[m], str):
            # Replace legacy method names with current ones
            method_str = variables[m].upper()
            method_str = legacy_map.get(method_str, method_str)

            try:
                variables[m] = class_dic[method][method_str]
            except KeyError:
                valid = list(class_dic[method].__members__.keys())
                raise ValueError(
                    f"Invalid {method}: {variables[m]!r}. "
                    f"Valid options: {valid}"
                )
    if len(variables) == 1:
        return variables[0]
    else:
        return variables


def validate_choice(value, choices, name):
    """Raise ValueError with a helpful option list when value is not in choices.

    Use for non-Enum string choices (e.g. well_type='gas'|'oil') where
    validate_methods does not apply.
    """
    if value not in choices:
        raise ValueError(
            f"Invalid {name}: {value!r}. Valid options: {list(choices)}"
        )
    return value