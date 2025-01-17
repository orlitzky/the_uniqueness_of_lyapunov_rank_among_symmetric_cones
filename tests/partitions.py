r"""
Functions for working with partitions.

In a few places in the paper, it suffices to consider only sums of
Lorentz cones, or maybe sums of Lorentz cones with exactly one
non-Lorentz factor. Sums of Lorentz cones are particularly nice
because they can be represented by partitions of the dimension you're
working in. For example, in dimension four, you can have ``L(4)``, or
``L(1)+L(3)``, or ``L(2)+L(2)``, or... you get the idea. Each possibility
is associated with a partition of ``4``, and vice-versa.

Eliminating isomorphic cones is also a bit easier when they are
represented by partitions. To eliminate any ambiguity arising from the
order of the factors, all one must do is sort the elements of the
corresponding partition; in fact, the :func:`partitions` function of
Kelleher and O'Sullivan does this already. There is only one more way
that ambiguity can arise, and that is from the equivalence of
``L(1)+L(1) == L(2)``. However this can be avoided quite easily by
ignoring any partitions that contain a ``2`` (which is fast, because
they are sorted).

For these reasons, we use partitions to test some of our
results. Here's a quick summary of this module's functions:

  * :func:`f` is the function from the paper that takes ``n`` and
    returns the Lyapunov rank of ``L(n)``.

  * :func:`partitions` is a slightly-modified version of the fast
    integer partitioning function of Kelleher and O'Sullivan.

  * :func:`partition_rank` computes the Lyapunov rank of the cone
    associated with a given partition.

  * :func:`partition_similacra` computes all similacra (represented as
    partitions) of the given cone (represented as a partition).

"""

from typing import Generator


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


def _partition_contains_two(p : list[int]) -> bool:
    r"""
    Return ``True`` if a _sorted_ partition (arising from
    the :func:`partitions` function) contains a ``2``.

    Parameters
    ----------

    p : list[int]
      The partition to check. Must be sorted least-to-greatest.

    Returns
    -------

    True if the given partition contains ``2``, and ``False`` otherwise.

    Examples
    --------

        >>> _partition_contains_two([1,1,1,1,3])
        False

        >>> _partition_contains_two([1,1,1,1,2,3])
        True

    If the partition isn't sorted, all bets are off::

        >>> _partition_contains_two([1,1,5,1,2])
        False

    """
    for z in p:
        if z == 2:
            return True
        if z >= 3:
            # Early return assuming the list is sorted
            return False
    return False


def partitions(n : int, entry_max : int | None = None, include_two : bool = True) -> Generator[list[int], None, None]:
    r"""
    Return all partitions of the integer ``n`` in ascending order.

    This is Jerome Kelleher's ``accel_asc`` function from
    https://jeromekelleher.net/category/combinatorics.html, modified
    in two ways:

      1. It takes an ``entry_max`` argument which limits the maximum size of
         any entry in a partition. Since the partitions are sorted, this is
         very easy to check before returning one of them.

      2. We can exclude partitions containing twos. One easy way to
         decrease the number of partitions under consideration is to
         normalize ``[2]`` to ``[1,1]``. Since a partition containing
         ``[1,1]`` in place of that ``[2]`` will be returned _anyway_,
         the easy way to do this is to drop all partitions containing
         a ``2``. We'll waste a tiny bit of space this way, but checking
         for ``2`` is a lot simpler than checking for ``1,1``.

    Parameters
    ----------

    n : int
      The integer to partition.

    entry_max : int | None
      An inclusive upper limit on the size of a partitions entries. If
      any entry in a partition exceeds this limit, it is
      omitted. Defaults to ``None`` (no limit).

    include_two : bool
      Whether or not to return partitions containing twos. (These are
      all equivalent to partitions containing ``1,1`` for the purposes
      of Lyapunov rank).

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
        >>> m = randint(1,20)
        >>> all( sum(p) == n for p in partitions(n,m) )
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

    If ``entry_max`` is one less than ``n``, that should only
    eliminate one partition (namely, ``[n]``)::

        >>> from random import randint
        >>> n = randint(1,20)
        >>> actual = len(tuple(partitions(n, n-1)))
        >>> expected = len(tuple(partitions(n))) - 1
        >>> actual == expected
        True

    Similarly, there's only one partition with all entries less than
    or equal to (i.e. equal to) one::

        >>> from random import randint
        >>> n = randint(1,20)
        >>> ps = tuple(partitions(n, 1))
        >>> len(ps)
        1
        >>> len(ps[0]) == n
        True

    We can exclude ``2`` from the results entirely as a cheap form of
    normalization::

        >>> list(partitions(2, include_two=False))
        [[1, 1]]

    This is a significant reduction in the number of partitions::

        >>> len(list(partitions(20)))
        627
        >>> len(list(partitions(20, include_two=False)))
        242

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
            if (entry_max is None) or a[k+1] <= entry_max:
                # "a" is sorted with the largest entries at the end
                if include_two or not _partition_contains_two(a[:k + 2]):
                    yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        if (entry_max is None) or a[k] <= entry_max:
            # "a" is sorted with the largest entries at the end
            if include_two or not _partition_contains_two(a[:k + 1]):
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


def partition_rank(p):
    r"""
    Return the Lyapunov rank of a sum of Lorentz factors whose
    dimensions are given by an integer partition.

    A direct sum of Lorentz cones is determined (almost) uniquely by
    the dimensions of its factors. The "almost" is because ``L(2)``
    and ``L(1) + L(1)`` are equal, but ``[1,1]`` and ``[2]`` are not.
    In any case, if we are given an integer partition that is intended
    to identify a direct sum of Lorentz cone (for example, computed by
    the :func:`partitions` function), then this function
    computes the Lyapunov rank of that direct sum.

    Parameters
    ----------

    p : [int]
      A partition of some integer, represented as a list
      of integers (whose sum if the one being partitioned).

    Returns
    -------

    An integer: the Lypaunov rank of the direct sum of Lorentz cones
    where the superscripts (i.e. the dimensions of the factors) are
    given by this partition.

    Examples
    --------

        >>> partition_rank([])
        0
        >>> partition_rank([0])
        0
        >>> partition_rank([1,2,3])
        7

    """
    return sum( map(f,p) )


def partition_similacra(p : list[int]) -> Generator[list[int], None, None]:
    r"""
    Find all partitions of the same integer as the given
    partition that have the same :func:`partitions.partition_rank` but
    are not equivalent.

    WARNING: the input partition should be sorted least to greatest,
    and should not contain ``2``. In other words, it should match the
    format returned by :func:`partitions`.

    We infer the integer from the given partition, and then compute
    the other partitions of it one-at-a-time. We skip partitions all
    of whose entries are less than the smallest entry in the given
    partition, since by Proposition 2, the Lyapunov rank associated
    with any such partition will be too small. We also skip partitions
    whose largest entry is equal to the smallest entry in the target
    partition, since the best we could hope for in that case is that
    all entries of both partitions are equal, and in that case the
    partitions are equal -- we want to toss those out regardless,
    because they are equivalent.

    After eliminating as many partitions as possible, we return
    the first with a rank that matches the given one.

    Parameters
    ----------

    p : list[int]
      The partition whose similacra you want. It should be SORTED,
      and should NOT CONTAIN ``2``.

    Returns
    -------

    Yields partitions of ``sum(p)`` that have the same
    :func:`partitions.partition_rank` as ``p``, but are not equivalent
    to it.

    Examples
    --------

    The nonnegative orthant will never have similacra::

        >>> list(partition_similacra([0]))
        []
        >>> list(partition_similacra([1]))
        []
        >>> list(partition_similacra([1,1,1,1,1]))
        []

    All other partitions of ``4`` have ranks that are too small (you
    can just try them all in your head)::

        >>> list(partition_similacra([4]))
        []

    The two similacra from Example 1::

        >>> next(partition_similacra([3,3,3,4]))
        [1, 1, 1, 1, 1, 1, 1, 1, 5]
        >>> sims = partition_similacra([1, 1, 1, 1, 1, 1, 1, 1, 5])
        >>> _ = next(sims); next(sims)
        [3, 3, 3, 4]

    If you ignore the warning and provide an unsorted partition or a
    partition containing ``2``, you may get wrong answers. Lack of
    sorting can eliminate valid similacra::

        >>> len(list(partition_similacra([3,3,3,4])))
        2
        >>> len(list(partition_similacra([4,3,3,3])))
        1

    And including ``2`` can produce equivalent partitions::

        >>> list(partition_similacra([2,5,13]))
        [[1, 1, 5, 13], [10, 10]]

    """
    target_rank = partition_rank(p)
    for q in partitions(sum(p), include_two=False):
        # The largest entry of q is its last, and the smallest entry
        # of p is its first.
        if q[-1] <= p[0]:
            continue
        if (partition_rank(q) == target_rank) and (not q == p):
            yield q
