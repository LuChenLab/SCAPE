from collections import defaultdict
from itertools import islice


# This two fucntion were used for wraper the nested dict not lambda function
def dict_list():
    return defaultdict(list)


def wrapdict():
    return defaultdict(dict_list)


def window(seq, n=2):
    """
    Return a sliding window over a list
    """
    it = iter(seq)
    result = list(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + list((elem,))
        yield result


class AttrDict(dict):

    def __init__(self):
        dict.__init__(self)

    def __setattr__(self, name, value):
        self[name] = value

    def __getattr__(self, name):
        return self[name]


class dotdict(dict):
    """
    dot.notation access to dictionary attributes
    """
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
