import pickle
import inspect

def load_vars(filename):
    caller_vars = inspect.stack()[1].frame.f_locals
    with open(filename, 'rb') as f:
        saved_vars = pickle.load(f)
    caller_vars.update(saved_vars)