r"""
Recursive functions to compute the admissible Lorentz ranks of a
given dimension, or to compute all symmetric cones (up to isomorphism)
in a given dimension.

These functions take both a ``db`` argument, and a ``sql`` toggle. If
you really want to affect the live database, you have to set ``sql``
to ``True``, and ``db`` to ``sql.LIVE_DATABASE``. These are not
default because we don't want to modify the live database by surprise
while testing (even though this shouldn't happen, ha ha).

If you run this module, on the other hand, it will begin to update the
live database::

    $ python compute.py
    computing dimension 99...

This however takes "forever" and will probably run your system out of
RAM. The bottleneck for both speed and space is putting all of the new
cones for a given dimension into a list and sorting them to eliminate
duplicates (specifically, serializations of isomorphic cones).
"""
import sqlite3

from cones import L,HR,HC,HH,HO
import signatures
import sql

def _admissible_lorentz_ranks(n : int, d : dict|None, db : str) -> tuple[int]:
    r"""
    Recursive implementation underlying :func:`admissible_lorentz_ranks`.

    This is necessary for that public function to have a nice user
    interface because the recursive bit will always pass a cache dict
    down to the next level, but we don't want users to have to pass in
    an empty dict to get started.

    If ``d`` is ``None`` instead of a dict, the SQL database ``db``
    will be consulted/updated instead.

    Parameters
    ----------

    n : int
      The dimension for which you'd like to know, what Lyapunov ranks
      are possible if we consider only Lorentz cone factors?

    d : dict|None
      Either a dict to cache the results in, or ``None`` if you want
      to use the SQL database ``db`` as a cache instead.

    db : str
      The name of the SQLite database to use as a cache (if ``d`` is
      ``None``).

    Returns
    -------

    tuple[int]
      A tuple of Lyapunov ranks that can be achieved in dimension
      ``n`` using only Lorentz cone factors, in no particular order.

    Examples
    --------

    The examples for :func:`admissible_lorents_ranks` demonstrate this
    indirectly, but we can check a few trivial cases by hand::

        >>> _admissible_lorentz_ranks(0, {}, "unused")
        (0,)
        >>> _admissible_lorentz_ranks(1, {}, "no database")
        (1,)

    This one will use a new, temporary database::

        >>> sql.new_database(sql.TEST_DATABASE)
        >>> sorted(_admissible_lorentz_ranks(3, None, sql.TEST_DATABASE))
        [3, 4]

    And the second time, it will be cached::

        >>> sorted(_admissible_lorentz_ranks(3, None, sql.TEST_DATABASE))
        [3, 4]

    To avoid surprises, ensure that the ``n = 0`` case gets added to the
    dictionary or database::

        >>> d = {}
        >>> _admissible_lorentz_ranks(2, d, None)
        (2,)
        >>> d
        {0: (0,), 1: (1,), 2: (2,)}

        >>> sql.new_database(sql.TEST_DATABASE)
        >>> _admissible_lorentz_ranks(2, None, sql.TEST_DATABASE)
        (2,)
        >>> sql.admissible_lorentz_ranks(0, sql.TEST_DATABASE)
        (0,)

    """
    if n == 1:
        # This function has side effects (updating the dictionary or
        # database) that are missed for n=0 because we never recurse
        # that far down. To avoid surprises, we call the theoretical
        # base case from the de facto base case to trigger its side
        # effects.
        _ = _admissible_lorentz_ranks(0, d, db)

    if d is None:
        if sql.have_lorentz_rank_dim(n, db=db):
            return sql.admissible_lorentz_ranks(n, db=db)
    else:
        if n in d:
            return d[n]

    s = set()
    for i in range(1,(n//2)+1):
        # We can stop at (n//2)+1 because afterwards, i passes n-i.
        # and the situation is symmetric.
        #
        # Compute s1 first (i.e. outside of the comprehension) because
        # we want to be sure that the lower values get computed before
        # we try to compute the higher ones.
        s1 = _admissible_lorentz_ranks(i, d, db)
        s2 = _admissible_lorentz_ranks(n-i, d, db)
        s.update( b1 + b2 for b1 in s1 for b2 in s2 )

    this_fn = ()
    if n != 2:
        # n=2 is the one case where f(n) = 1 + 1 + ... + 1 (n times)
        # and the cone is reducible, so f(n) would wind up in the list
        # twice, once for L(n) and once for L(1) + L(1).
        this_fn = (signatures.f(n),)
    result = tuple(s) + this_fn

    if d is None:
        sql.insert_lorentz_ranks(n, result, db=db)
    else:
        d[n] = result
    return result


def admissible_lorentz_ranks(n : int, sql : bool = False, db : str = sql.TEST_DATABASE) -> tuple[int]:
    r"""
    Compute admissible Lyapunov ranks for sums of Lorentz cones
    of total dimension ``n``.

    This operates recursively and caches the results at each step (to
    make future steps faster). It can also use a SQL database (instead
    of the default python dict) as a cache.

    Parameters
    ----------

    n : int
      The dimension for which you'd like to know, what Lyapunov ranks
      are possible if we consider only Lorentz cone factors?

    sql : bool, default=False
      Whether or not to use the SQL database, or start fresh.

    db : str, default=TEST_DATABASE
      The name of the SQLite database to use (if ``sql`` is ``True``).

    Returns
    -------

    tuple[int]
      The "set" of Lyapunov ranks that can be achieved in dimension
      ``n`` using only Lorentz cone factors. A tuple is used instead
      of a set because msgpack knows what to do with a tuple, unlike
      a set.

    Examples
    --------

    Low-dimensional examples::

        >>> admissible_lorentz_ranks(0)
        (0,)
        >>> admissible_lorentz_ranks(1)
        (1,)
        >>> admissible_lorentz_ranks(2)
        (2,)
        >>> sorted(admissible_lorentz_ranks(3))
        [3, 4]

    The SQL values should agree with the ones we compute. We can
    check this in a reasonable amount of time up to ``n = 75``. We
    sort the results before comparing them because, despite our use of
    tuples, the order that they wind up in is not meaningful::

        >>> def check(n):
        ...     if not sql.have_lorentz_rank_dim(n):
        ...         # don't fail if we're just missing the data
        ...         return True
        ...     actual = sorted(sql.admissible_lorentz_ranks(n))
        ...     expected = sorted(admissible_lorentz_ranks(n))
        ...     return (actual == expected)
        >>> check(0)
        True
        >>> check(1)
        True
        >>> check(2)
        True
        >>> check(3)
        True
        >>> from random import randint
        >>> n = randint(1,76)
        >>> check(n)
        True

    Since we don't recurse all the way down to ``n = 0``, some care is
    needed to make sure that we don't assume the existence of that row
    based on the presence of rows for larger dimensions::

        >>> sql.new_database(sql.TEST_DATABASE)
        >>> sorted(admissible_lorentz_ranks(4, True))
        [4, 5, 7]
        >>> admissible_lorentz_ranks(0, True)
        (0,)

    """
    # The implementation of this function _always_ uses a cache,
    # the only question is, whether or not the cache will be
    # a python dict that gets passed around, or an implicit
    # SQL database.
    d = {}
    if sql:
        d = None

    # Now just run the real, recursive implementation.
    return _admissible_lorentz_ranks(n, d, db)


def _irreducible_cones_of_dim(n : int) -> tuple:
    r"""
    Return a tuple of irreducible cones in dimension ``n``.

    Outside of dimension two, there is always a Lorentz cone in
    dimension ``n``, but there may not be any others. This is "up to
    isomorphism," so, for example, you won't get ``HR(2)`` in
    dimension three.

    Parameters
    ----------

    n : int
      The dimension in which you'd like the irreducible cones.

    Returns
    -------

    A tuple consisting of all irreducible cones in dimension ``n``.

    Examples
    --------

    There are no irreducible cones of dimension two::

        >>> _irreducible_cones_of_dim(0)
        (L(0),)
        >>> _irreducible_cones_of_dim(1)
        (L(1),)
        >>> _irreducible_cones_of_dim(2)
        ()

    Isomorphic results are not returned::

        >>> L(3).dim == HR(2).dim == 3
        True
        >>> _irreducible_cones_of_dim(3)
        (L(3),)

    Larger factors appear where we think they will::

        >>> HR(3) in _irreducible_cones_of_dim(6)
        True
        >>> HO(3) in _irreducible_cones_of_dim(27)
        True

    """
    s = []
    if n != 2:
        s.append(L(n))

    if n >= 27:    # HO(3).dim
        c = HO.in_dim(n)
        if c:
            s.append(c)
    elif n >= 15:  # HH(3).dim
        c = HH.in_dim(n)
        if c:
            s.append(c)
    elif n >= 9:   # HC(3).dim
        c = HC.in_dim(n)
        if c:
            s.append(c)
    elif n >= 6:   # HR(3).dim
        c = HR.in_dim(n)
        if c:
            s.append(c)

    return tuple(s)


def _merge_factors(a : int|tuple[int], b : int|tuple[int]) -> tuple[int]:
    r"""
    Merged two serialized cones into a third.

    We can do this efficiently because the sort order for cone
    factors is already based on their serializations. As a result,
    we don't need to deserialize and reserialize; we can just
    combine and sort the serialized representations directly.

    Parameters
    ----------

    a : int|tuple[int]
      The first cone, serialized (as either an int or a tuple of
      ints).
    b : int|tuple[int]
      The second cone, serialized (as either an int or a tuple of
      ints).

    Returns
    -------

    A tuple of ints representing a direct sum. If the two inputs ``a``
    and ``b`` were deserialized and then combined into a
    :class:`DirectSum`, the serialization of that direct sum is what
    we return.

    Examples
    --------

    A simple example with two irreducible factors::

        >>> from cones import SymmetricCone
        >>> a = L(4).serialize()
        >>> b = HR(3).serialize()
        >>> c = _merge_factors(a, b); c
        (32, 41)
        >>> SymmetricCone.deserialize(c)
        HR(3) + L(4)

    A random example showing that this does what we think it does,
    i.e. circumvents deserialization/reserialization accurately::

        >>> from cones import DirectSum, SymmetricCone, random_cone
        >>> K = random_cone()
        >>> a = K.serialize()
        >>> J = random_cone()
        >>> b = J.serialize()
        >>> expected = DirectSum([K,J])
        >>> actual = SymmetricCone.deserialize(_merge_factors(a,b))
        >>> actual == expected
        True

    """
    # computing/comparing the type to int is actually
    # a bit faster on average than isinstance
    if type(a) == int:
        a = (a,)
    if type(b) == int:
        b = (b,)

    # We have to re-sort the factors, because that's what the direct
    # sum constructor would do to ensure that no duplicates sneak
    # in. The serializations (s1,s2) and (s2,s1) represent the same
    # cone! The caller is responsible for deduplicating them.
    return tuple(sorted(a+b))


def _dim_ranks_cones(n : int, d : dict|None, db : str, progress : bool) -> dict:
    r"""
    Recursive implementation underlying :func:`dim_ranks_cones`.

    This is necessary for that public function to have a nice user
    interface because the recursive bit will always pass a cache dict
    down to the next level, but we don't want users to have to pass in
    an empty dict to get started.

    If ``d`` is ``None`` instead of a dict, the SQL database ``db``
    will be consulted/updated instead.

    Parameters
    ----------

    n : int
      The dimension for which you want the rank => cones map.

    d : dict|None
      Either a dict to cache the results in, or ``None`` if you want
      to use the SQL database ``db`` as a cache instead.

    db : str
      The name of the SQLite database to use as a cache (if ``d`` is
      ``None``).

    progress : bool
      Whether or not to print a dot "." as the function progresses.
      This isn't tested or anything, but it's real nice to see things
      moving when computing the (n+1)st dimension takes two days.

    Returns
    -------

    dict
      A dictionary whose keys are all possible Lyapunov ranks in
      dimension ``n``, and whose values are tuples of (serialized)
      cones in dimension ``n`` having that Lyapunov rank.

    Examples
    --------

    The examples for :func:`dim_ranks_cones` demonstrate this
    indirectly, but we can check a few trivial cases by hand::

        >>> _dim_ranks_cones(0, {}, "unused", False)
        {0: (1,)}
        >>> _dim_ranks_cones(1, {}, "no database", False)
        {1: (11,)}

    This one will use a new, temporary database::

        >>> sql.new_database(sql.TEST_DATABASE)
        >>> _dim_ranks_cones(3, None, sql.TEST_DATABASE, False)
        {3: ((11, 11, 11),), 4: (31,)}

    And the second time, it will be cached::

        >>> _dim_ranks_cones(3, None, sql.TEST_DATABASE, False)
        {3: ((11, 11, 11),), 4: (31,)}

    To avoid surprises, ensure that the ``n = 0`` case gets added to the
    dictionary or database::

        >>> d = {}
        >>> _ = _dim_ranks_cones(2, d, None, False)
        >>> d[0]
        {0: (1,)}

        >>> sql.new_database(sql.TEST_DATABASE)
        >>> _ = _dim_ranks_cones(2, None, sql.TEST_DATABASE, False)
        >>> sql.dim_ranks_cones(0, sql.TEST_DATABASE)
        {0: (1,)}

    """
    if n == 1:
        # This function has side effects (updating the dictionary or
        # database) that are missed for n=0 because we never recurse
        # that far down. To avoid surprises, we call the theoretical
        # base case from the de facto base case to trigger its side
        # effects.
        _ = _dim_ranks_cones(0, d, db, progress)

    if d is None:
        if sql.have_cone_dim(n, db=db):
            return sql.dim_ranks_cones(n, db=db)
    else:
        if n in d:
            return d[n]

    # The dict for this n. It will either be inserted as d[n],
    # or put into the SQL database instead.
    d_n = {}

    # We partition "n" ourselves here. Basically, we split n into (i,
    # n-i), and then recurse into each of them. Every partition of "n"
    # arises from a partition of "i" plus a partition of "n-i". We can
    # stop at n//2 here because the situation is symmetric: once i
    # passes n-i, the resulting set is going to be the same; (5,2)
    # gives the same result as (2,5).
    for i in range(1,(n//2)+1):
        s1 = _dim_ranks_cones(i, d, db, progress)
        if progress:
            print(".", end="", flush=True)
        s2 = _dim_ranks_cones(n-i, d, db, progress)
        if progress:
            print(".", end="", flush=True)

        # the ranks possible in dim=n are the sums of ranks possible
        # in dim=i and dim=(n-i)
        for r1 in s1:
            for r2 in s2:
                s = set( _merge_factors(b1,b2)
                         for b1 in s1[r1]
                         for b2 in s2[r2] )
                r = r1 + r2
                if r in d_n:
                    # there's more than one way to add up to r!
                    d_n[r].update(s)
                else:
                    d_n[r] = s

    for c in _irreducible_cones_of_dim(n):
        if c.rank in d_n:
            d_n[c.rank].add(c.serialize())
        else:
            d_n[c.rank] = {c.serialize()}

    for r in d_n:
        # for serialization we want tuples, not sets
        d_n[r] = tuple(d_n[r])

    if d is None:
        sql.insert_cones(n, d_n, db=db)
    else:
        d[n] = d_n
    return d_n


def dim_ranks_cones(n : int, sql : bool = False, db : str = sql.TEST_DATABASE, progress : bool = False) -> dict:
    r"""
    Compute a rank => cones map for all cones of dimension ``n``.

    Parameters
    ----------
    n : int
      The dimension for which you want the rank => cones dict.

    sql : bool, default=False
      Whether or not to use the SQL database, or start fresh.

    db : str, default=TEST_DATABASE
      The name of the SQLite database to use (if ``sql`` is ``True``).

    progress : bool, default=False
      Whether or not to print a dot "." as the function progresses.
      This isn't tested or anything, but it's real nice to see things
      moving when computing the (n+1)st dimension takes two days.

    Examples
    --------

    Easy low-dimensional cases::

        >>> dim_ranks_cones(3)
        {3: ((11, 11, 11),), 4: (31,)}

        >>> dim_ranks_cones(4)
        {4: ((11, 11, 11, 11),), 5: ((11, 31),), 7: (41,)}

        >>> dim_ranks_cones(5)
        {5: ((11, 11, 11, 11, 11),), 6: ((11, 11, 31),), 8: ((11, 41),), 11: (51,)}

    The precomputed values should agree with the ones we compute::

        >>> def check(n):
        ...     d1 = sql.dim_ranks_cones(n)
        ...     d2 = dim_ranks_cones(n)
        ...     return ( all( set(d1[r]) == set(d2[r]) for r in d1 )
        ...              and sorted(d1.keys()) == sorted(d2.keys()) )
        >>> from random import randint
        >>> n = randint(0,40)
        >>> check(n)
        True

    Since we don't recurse all the way down to ``n = 0``, some care is
    needed to make sure that we don't assume the existence of that row
    based on the presence of rows for larger dimensions::

        >>> sql.new_database(sql.TEST_DATABASE)
        >>> _ = dim_ranks_cones(4, True)
        >>> dim_ranks_cones(0, True)
        {0: (1,)}

    """
    # The implementation of this function _always_ uses a cache, the
    # only question is, whether or not the cache will be a python dict
    # that gets passed around, or an implicit SQL database.
    d = {}
    if sql:
        d = None
    # Now just run the real, recursive implementation.
    return _dim_ranks_cones(n, d, db, progress)



if __name__ == "__main__":
    # if executed, we start computing more cones
    mcd = sql.max_cone_dim()
    if mcd:
        n = mcd + 1
    else:
        # there won't be a maximum dimension if the DB is empty
        n = 0

    while True:
        print(f"computing dimension {n}", end="", flush=True)
        _ = dim_ranks_cones(n, True, sql.LIVE_DATABASE, True)
        print(" done.")
        n += 1
