r"""
Compute possible signatures for symmetric cones.
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


def _partitions(n : int) -> list[list[int]]:
    r"""
    Return all partitions of the integer ``n`` in ascending order.

    This is Jerome Kelleher's ``accel_asc`` function from
    https://jeromekelleher.net/category/combinatorics.html.

    Parameters
    ----------

    n : int
      The integer to partition.

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
    >>> all( sum(p) == n for p in _partitions(n) )
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
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
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

    Some comparisons with the recursive algorithm::

    >>> all( _direct_lorentz_ranks(k) == admissible_lorentz_ranks(k)
    ...      for k in [7,12,15,19,23,30] )
    True

    """
    # dedupe
    return tuple(set( (sum(f(p_k) for p_k in p) )
                      for p in _partitions(n) ))
