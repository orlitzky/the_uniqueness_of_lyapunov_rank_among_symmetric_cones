r"""
Convert the all_cones_of_dimension() cache to one for
dim_ranks_cones().

This is nice to have as a separate program because the
all_cones_of_dimension() computation is essentially sequential,
whereas deserialization can happen in parallel (yes, instances are
cached, but there's not a strict dependence of the bigger cones on the
smaller ones).
"""
import gzip, msgpack
from multiprocessing import Lock
from multiprocessing.pool import ThreadPool
from cones import *

def compute(n, _acod_cache, _drc_cache, l):
    print(f"computing _drc_cache[{n}] from _acod_cache[{n}]... ",
          flush=True)

    # work outside of the real dict to avoid locking
    # as long as possible
    _drc_cache_n = {}
    for c_ser in _acod_cache[n]:
        # get the lock because deserialization affects instance
        # cache dicts
        l.acquire()
        c = SymmetricCone.deserialize(c_ser)
        l.release()
        if c.rank in _drc_cache_n:
            _drc_cache_n[c.rank].append(c_ser)
        else:
            _drc_cache_n[c.rank] = [c_ser]
    l.acquire()
    _drc_cache[n] = _drc_cache_n
    l.release()

    print(f"done {n}", flush=True)

if __name__ == "__main__":
    with gzip.open("_acod_cache.pck.gz", "rb") as f:
        _acod_cache = msgpack.unpack(f,
                                     strict_map_key=False,
                                     use_list=False)


    # update an existing _drc_cache if we can
    try:
        with gzip.open("_drc_cache.pck.gz", "rb") as f:
            _drc_cache = msgpack.unpack(f,
                                         strict_map_key=False,
                                         use_list=False)
    except FileNotFoundError:
        # The cache dict that we'll eventually write to disk.
        _drc_cache = {
            0: { 0 : () },
            1: { 1 : (11,) },
            2: { 2: ((11,11),) }
        }

    n_min = max(_drc_cache.keys()) + 1
    n_max = max(_acod_cache.keys())

    # passed to each thread so that we don't write to
    # any dictionaries simultaneously
    l = Lock()

    with ThreadPool(40) as p:
        results = [ p.apply_async(compute, (n, _acod_cache, _drc_cache, l))
                    for n in range(n_min, n_max+1) ]
        _ = [ r.get() for r in results ]

    print("serializing...", end="", flush=True)
    with gzip.open("_drc_cache.pck.gz", "wb") as f:
        msgpack.pack(_drc_cache, f)
    print("done")
