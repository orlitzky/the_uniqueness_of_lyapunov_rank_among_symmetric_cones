r"""
Partition-related stuff.
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


from typing import Iterator
def _partition_contains_two(p : Iterator[int]) -> bool:
    r"""
    Return ``True`` if a _sorted_ partition (arising from
    the :func:`partitions` function) contains a ``2``.

    Parameters
    ----------

    p : Iterator[int]
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


def partitions(n : int, entry_max : int = None, include_two : bool = True) -> list[list[int]]:
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

    entry_max : int
      An inclusive upper limit on the size of a partitions entries. If
      any entry in a partition exceeds this limit, it is omitted.

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


def partitions_equivalent(p,q):
    r"""
    Determine whether two partitions represent the same
    sum-of-Lorentz cones.

    The order of factors in a direct sum does not matter, so one thing
    this function does is sort the elements of each partition. But we
    also must take care that ``RN(2)`` and ``L(2)`` are the same cone.

    Parameters
    ----------

    p,q : [int]
      Integer partitions representing sums of Lorentz cones.

    Returns
    -------

    ``True`` if ``p`` and ``q`` represent the same sum-of-Lorentz
    cone, and ``False`` otherwise.

    Examples
    --------

    Small batch, hand-crafted examples::

        >>> partitions_equivalent([], [0])
        True
        >>> partitions_equivalent([1], [2])
        False
        >>> partitions_equivalent([0,0], [])
        True
        >>> partitions_equivalent([1,2,3], [1,2,2,1])
        False
        >>> partitions_equivalent([1,1,2], [1,1,1,1])
        True
        >>> partitions_equivalent([2,1,2], [1,1,1,1,1])
        True
        >>> partitions_equivalent([1,1,3], [2,3])
        True
        >>> partitions_equivalent([1,1,3], [3,2])
        True
        >>> partitions_equivalent([1,3,3], [2,2,3])
        False

    Being equivalent is a symmetric relationship::

        >>> from random import randint
        >>> n = randint(1,10)
        >>> ps = list(partitions(n))
        >>> all(
        ...   partitions_equivalent(p,q)
        ...   ==
        ...   partitions_equivalent(q,p)
        ...   for p in ps
        ...   for q in ps
        ... )
        True

    """
    if (sum(p) != sum(q)):
        return False

    # sort and remove zeros (which shouldn't be there in the first
    # place)
    p = sorted(i for i in p if not i == 0)
    q = sorted(j for j in q if not j == 0)

    # Remove all factors of size >= 3 in p from both p and q. Use
    # indices instead of a "for foo in bar" loop because we're going
    # to be deleting items from the list we're iterating over.
    idx = 0
    len_p = len(p)
    while idx < len_p:
        if p[idx] <= 2:
            # skip it
            idx += 1
        elif p[idx] not in q:
            # Since p[idx] >= 3, if it's missing from q, they're not
            # isomorphic. (Without p[idx] >= 3 this doesn't work,
            # because for example [1,1] and [2] have no elements in
            # common.)
            return False
        else:
            # Otherwise, remove this element from both p and q, but
            # don't increment idx, because they'll all shift down by
            # one.
            q.remove(p[idx])
            del(p[idx])
            len_p -= 1

    # Now what's left in p is its 1,2 elements; and what's left in q
    # is whatever 1,2 elements it had plus any elements >= 3 that
    # were not in p. We can repeat in the opposite direction.
    idx = 0
    len_q = len(q)
    while idx < len_q:
        if q[idx] <= 2:
            # skip it
            idx += 1
        elif q[idx] not in p:
            return False
        else:
            p.remove(q[idx])
            del(q[idx])
            len_q -= 1

    # Now both p and q should have only 1s and 2s in them. They should
    # still sum to the same value if they are isomorphic.
    if (not p) or (not q):
        # one of them's empty, they both had better be
        return (p == q)

    # Both nonempty, this is safe
    if max(p) > 2 or max(q) > 2:
        raise ValueError("elements of size >= 2 left in partition")

    # Otherwise, so long as they still add up to the same size, the
    # 1-dim and 2-dim factors can all be grouped.
    return sum(p) == sum(q)
