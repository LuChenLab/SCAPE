import random
import numpy as np

def replace_with_closest(ref_arr, query_arr):
    n = len(query_arr)
    res = np.zeros(n)
    for i in range(n):
        tmp_ind = np.argmin(abs(ref_arr - query_arr[i]))
        res[i] = ref_arr[tmp_ind]
    return res


def gen_k_arr(K, n):
    """
    Arguments:
        K {int} -- [apa numbers]
        n {int} -- [trial numbers]
    """

    def random_sel(K, trial=200):
        count_index = 0
        pool = np.arange(K)
        last = None
        while count_index < trial:
            count_index += 1
            random.shuffle(pool)
            if pool[0] == last:
                swap_with = random.randrange(1, len(pool))
                pool[0], pool[swap_with] = pool[swap_with], pool[0]
            for item in pool:
                yield item

            last = pool[-1]

    if K <= 1:
        return np.repeat(K - 1, n)
    else:
        k_lst = list(random_sel(K, trial=n))
        return np.array(k_lst)


def rle(x):
    """
    https://gist.github.com/nvictus/66627b580c13068589957d6ab0919e66#file-runlength-py-L13
    """
    where = np.flatnonzero
    x = np.asarray(x)
    n = len(x)
    if n == 0:
        return (np.array([], dtype=int),
                np.array([], dtype=int),
                np.array([], dtype=x.dtype))

    starts = np.r_[0, where(~np.isclose(x[1:], x[:-1], equal_nan=True)) + 1]
    lengths = np.diff(np.r_[starts, n])
    values = x[starts]

    return starts, lengths, values


# For time consulting testing
def func_line_time(f):
    @wraps(f)
    def decorator(*args, **kwargs):
        func_return = f(*args, **kwargs)
        lp = LineProfiler()
        lp_wrap = lp(f)
        lp_wrap(*args, **kwargs) 
        lp.print_stats() 
        return func_return 
    return decorator 
