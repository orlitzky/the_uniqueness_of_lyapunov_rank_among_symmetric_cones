r"""
Compute possible signatures for symmetric cones.

This module is mostly obsolete, but it can still be used as a sanity
check for the recursive/cached implementations of the same procedure.
"""

def f(n : int) -> int:
    r"""
    Compure the Lyapunov rank of the Lorentz cone in ``n``
    dimensions. This function was imaginatively called ``f`` in
    Proposition 8 and Theorem 4.

    Parameters
    ----------

    n : int
      The dimension of the Lorentz cone whose Lyapunov rank you
      want. Should be nonnegative unless you want a nonsense
      answer.

    Returns
    -------

    int
      The Lyapunov rank of the Lorentz cone in ``n`` dimensions

    Examples
    --------

    Small ``n`` are easy to check by hand::

        >>> f(0)
        0
        >>> f(1)
        1
        >>> f(2)
        2
        >>> f(3)
        4

    """
    if n == 0:
        return 0
    return (n**2 - n + 2) // 2


def partitions(n : int, entry_max : int = None) -> list[list[int]]:
    r"""
    Return all partitions of the integer ``n`` in ascending order.

    This is Jerome Kelleher's ``accel_asc`` function from
    https://jeromekelleher.net/category/combinatorics.html, modified
    to take an ``entry_max`` argument which limits the maximum size of
    any entry in a partition.

    Parameters
    ----------

    n : int
      The integer to partition.

    entry_max : int
      An inclusive upper limit on the size of a partitions entries. If
      any entry in a partition exceeds this limit, it is omitted.

    Returns
    -------

    list[list[int]]
      A list of partitions of ``n``. Each "partition" it
      itself a list of integers that sums to ``n``.

    Examples
    --------

    Each "partition" should actually partition ``n``::

        >>> from random import randint
        >>> n = randint(1,20)
        >>> all( sum(p) == n for p in partitions(n) )
        True

    If ``entry_max`` is larger than the integer we're partitioning, it
    should have no effect::

        >>> from random import randint
        >>> n = randint(1,15)
        >>> tuple(partitions(n)) == tuple(partitions(n,25))
        True

    Check that the ``entry_max`` is respected::

        >>> from random import randint
        >>> n = randint(1,15)
        >>> entry_max = randint(0,n)
        >>> all( z <= entry_max
        ...      for p in partitions(n, entry_max)
        ...      for z in p )
        True

    """
    a = [0 for i in range(n + 1)]
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            if (entry_max is None) or all(z <= entry_max for z in a[:k + 2]):
                yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        if (entry_max is None) or all(z <= entry_max for z in a[:k + 1]):
            yield a[:k + 1]


def _direct_lorentz_ranks(n : int) -> tuple:
    r"""
    Generate admissible Lyapunov ranks for all sums of Lorentz
    cones whose dimensions add up to ``n``.

    This is a direct (and much slower) computation using a different
    strategy than :func:`admissible_lorentz_ranks`.

    Parameters
    ----------

    n : int
      The dimension to compute signatures for.

    Returns
    -------

    tuple[int]
      A tuple of Lypaunov ranks.

    Examples
    --------

    In dimensions two and fewer, the Lorentz cone is the nonnegative
    orthant::

        >>> _direct_lorentz_ranks(0)
        (0,)
        >>> _direct_lorentz_ranks(1)
        (1,)
        >>> _direct_lorentz_ranks(2)
        (2,)

    Some comparisons with the recursive algorithm / database::

        >>> import compute
        >>> import sql
        >>> all( _direct_lorentz_ranks(k)
        ...      ==
        ...      sql.admissible_lorentz_ranks(k)
        ...      ==
        ...      compute.admissible_lorentz_ranks(k)
        ...      for k in [7,12,15,19,23] )
        True

    """
    # go from list -> set -> tuple to deduplicate them
    return tuple(set( (sum(f(p_k) for p_k in p) )
                      for p in partitions(n) ))
