"""Convert pyrestoolbox return values to JSON-serializable structures."""

import dataclasses
import math
from enum import Enum


def to_jsonable(obj):
    """Recursively convert obj to JSON-serializable Python primitives.

    numpy scalars/arrays, pandas DataFrames/Series, dataclasses, Enums and
    nested containers are handled. Non-finite floats become None (JSON has
    no NaN/inf).
    """
    if obj is None or isinstance(obj, (bool, str)):
        return obj
    if isinstance(obj, int):
        return int(obj)
    if isinstance(obj, float):  # np.float64 subclasses float - cast to plain float
        return float(obj) if math.isfinite(obj) else None
    if isinstance(obj, Enum):
        return obj.name
    if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
        return {f.name: to_jsonable(getattr(obj, f.name)) for f in dataclasses.fields(obj)}
    # numpy
    try:
        import numpy as np
        if isinstance(obj, np.generic):
            return to_jsonable(obj.item())
        if isinstance(obj, np.ndarray):
            return [to_jsonable(v) for v in obj.tolist()]
    except ImportError:
        pass
    # pandas
    try:
        import pandas as pd
        if isinstance(obj, pd.DataFrame):
            return {
                "columns": [str(c) for c in obj.columns],
                "records": to_jsonable(obj.to_dict(orient="records")),
            }
        if isinstance(obj, pd.Series):
            return to_jsonable(obj.tolist())
    except ImportError:
        pass
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [to_jsonable(v) for v in obj]
    return str(obj)
