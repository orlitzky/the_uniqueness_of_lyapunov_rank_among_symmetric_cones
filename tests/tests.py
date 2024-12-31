r"""
Exhaustive integer partition calculations to confirm results.

Test the expressions derived in Proposition 1::

    >>> from sympy import expand, floor, simplify, symbols
    >>> m,n,k = symbols("m,n,k", integer=True, positive=True)

Sympy doesn't know (for example) that n^2 + n is even, so we have to
fix the "floor" that it inserts everywhere we use integer
division-by-two on an even expression::

    >>> K = L(n)
    >>> J1 = L(k)
    >>> J2 = L(n-k)
    >>> J = DirectSum([J1,J2], False)  # can't sort symbolics
    >>> actual = fix_floor(K.rank - J.rank)
    >>> expected = (n-k)*k - 1
    >>> simplify(actual - expected)
    0

    >>> K = L(27)
    >>> J = HO(3)
    >>> K.dim == J.dim
    True
    >>> K.rank > J.rank
    True

    >>> K = L((m**2 + m)/2)
    >>> J = HR(m)
    >>> fix_floor(K.dim - J.dim)
    0
    >>> expand(fix_floor(K.rank - J.rank))
    m**4/8 + m**3/4 - 9*m**2/8 - m/4 + 1

    >>> K = L(m**2)
    >>> J = HC(m)
    >>> fix_floor(K.dim - J.dim)
    0
    >>> fix_floor(K.rank - J.rank)
    m**4/2 - 5*m**2/2 + 2

    >>> K = L(2*m**2 - m)
    >>> J = HH(m)
    >>> fix_floor(K.dim - J.dim)
    0
    >>> expand(fix_floor(K.rank - J.rank))
    2*m**4 - 2*m**3 - 9*m**2/2 + m/2 + 1

We verify Proposition 8: if ``n > 30``, we never get similacra. The
precomputed signatures are used for this so that it completes in a
reasonable time::

    >>> from sql import admissible_lorentz_ranks
    >>> nmax = 241
    >>> nmin = 31
    >>>
    >>> all(
    ...   (17+L(n).rank)
    ...   not in admissible_lorentz_ranks(9+n)
    ...   for n in range(nmin,nmax+1)
    ... )
    True

Confirm the table for ``n <= 30``::

    >>> d = {
    ...   2:  [5,3,3],
    ...   3:  [4,4,4],
    ...   4:  [6,3,1,1,1,1],
    ...   5:  [6,4,3,1],
    ...   6:  [7,4,1,1,1,1],
    ...   7:  [8,3,3,1,1],
    ...   8:  [9,3,1,1,1,1,1],
    ...   9:  [10,1,1,1,1,1,1,1,1],
    ...   10: [9,7,3],
    ...   15: [14,8,1,1],
    ...   18: [14,13],
    ...   21: [19,11],
    ...   22: [21,9,1],
    ...   30: [29,10] }
    >>>
    >>> def check(n):
    ...     J = DirectSum([ HC(3), L(n) ])
    ...     if n in d:
    ...         K = DirectSum(tuple( L(i) for i in d[n] ))
    ...         return(
    ...           J.rank in admissible_lorentz_ranks(J.dim)
    ...           and
    ...           K.signature() == J.signature()
    ...         )
    ...     else:
    ...         return J.rank not in admissible_lorentz_ranks(J.dim)
    >>>
    >>> all( check(n) for n in range(31) )
    True

To be extra sure, we check the table for Proposition 8 using the
similacra method as well::

    >>> n_with_similacra = []
    >>>
    >>> for n in range(31):
    ...     K = DirectSum([HC(3),L(n)])
    ...     if len(K.similacra()) != 0:
    ...         n_with_similacra.append(n)
    >>>
    >>> n_with_similacra
    [2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 18, 21, 22, 30]

Check the symbolic identity in Lemma 6 for the derivative of the
g-delta function::

    >>> from sympy import diff, symbols
    >>> x,d = symbols("x,d", integer=True, positive=True)
    >>> g = fix_floor(L(x).rank - L(x-d).rank)
    >>> diff(g, x)
    d

Check the symbolic identity in Lemma 6 for the Lyapunov rank of
``L(n-1)`` in terms of that of ``L(n)``::

    >>> from sympy import simplify, symbols
    >>> n = symbols("n", integer=True, positive=True)
    >>> lhs = L(n-1).rank
    >>> rhs = L(n).rank - (n-1)
    >>> simplify(fix_floor(lhs - rhs))
    0

Check implication (3) in Lemma 6::

    >>> f = lambda x: L(x).rank
    >>> all( f(n-1) + f(1+dimK) >= f(n-1-d) + f(1+dimK+d)
    ...      for d in range(100)
    ...      for dimK in range(100)
    ...      for n in range(2 + dimK + d, 100) )
    True

Check the relationships between the lower bounds on ``n``. This isn't
stated anywhere in the paper, but it shows that no bound is implied by
the other (each bound can be strictly the largest)::

    >>> K = L(2)
    >>> lowerbound1(K) < lowerbound2(K) < lowerbound3b(K)
    True

    >>> K = RN(5)
    >>> lowerbound1(K) < lowerbound3b(K) < lowerbound2(K)
    True

    >>> K = L(5)
    >>> lowerbound3b(K) < lowerbound2(K) < lowerbound1(K)
    True

Verify the cases mentioned explicitly in Proposition 11. First, the
``m != 2`` cases where there are no similacra::

    >>> m = 4
    >>> K = L(m)
    >>> n = 5
    >>> (n >= lowerbound1(K), n >= lowerbound2(K), n >= lowerbound3b(K))
    (True, False, True)
    >>> DirectSum([K,L(n)]).similacra()
    ()

    >>> m = 3
    >>> K = L(m)
    >>> n = 4
    >>> (n >= lowerbound1(K), n >= lowerbound2(K), n >= lowerbound3b(K))
    (True, False, False)
    >>> DirectSum([K,L(n)]).similacra()
    ()
    >>> n = 3
    >>> (n >= lowerbound1(K), n >= lowerbound2(K), n >= lowerbound3b(K))
    (True, False, False)
    >>> DirectSum([K,L(n)]).similacra()
    ()

    >>> m = 1
    >>> K = L(m)
    >>> n = 2
    >>> (n >= lowerbound1(K), n >= lowerbound2(K), n >= lowerbound3b(K))
    (True, False, False)
    >>> DirectSum([K,L(n)]).similacra()
    ()
    >>> n = 3
    >>> (n >= lowerbound1(K), n >= lowerbound2(K), n >= lowerbound3b(K))
    (True, True, False)
    >>> DirectSum([K,L(n)]).similacra()
    ()
    >>> n = 4
    >>> (n >= lowerbound1(K), n >= lowerbound2(K), n >= lowerbound3b(K))
    (True, True, False)
    >>> DirectSum([K,L(n)]).similacra()
    ()

And now the ``m == 2`` case where there is exactly one counterexample::

    >>> m = 2
    >>> K = L(m)
    >>> lowerbound1(K)
    2
    >>> DirectSum([K,L(2)]).similacra()
    ()
    >>> DirectSum([K,L(3)]).similacra()
    ()
    >>> DirectSum([K,L(4)]).similacra()
    (HR(3),)

Check Proposition 11 directly::

    >>> def check(m,n):
    ...     K = DirectSum([L(m),L(n)])
    ...     return not K.similacra()
    >>> from sql import max_cone_dim
    >>> # need m+n <= max_cone_dim()
    >>> m_max = max_cone_dim() // 2
    >>> all( check(m,n)
    ...      for m in range(1, m_max)
    ...      for n in range(lowerbound1(L(m)), max_cone_dim() - m)
    ...      if m != 2 )
    True

Check Giovanni's Theorem 4, as far as we can::

    >>> from cones import DirectSum, L
    >>> from sql import max_cone_dim
    >>> n_max = max_cone_dim() // 2
    >>> n_min = 1
    >>> n_without_similacra = []
    >>> for n in range(n_min, n_max+1):
    ...     K = DirectSum([L(n)]*2)
    ...     if not K.similacra():
    ...         n_without_similacra.append(n)
    >>> n_without_similacra
    [1, 2, 3, 5, 6, 7, 11, 12, 13, 18]


In Theorem 4, we "replace" the non-Lorentz irreducible factors with
sums of Lorentz cones. Here we confirm that those sums have the
correct dimensions, and Lyapunov ranks that dominates the Lyapunov
ranks of the things they replace. The first example we give is for the
complex PSD cones::

    >>> from sympy import symbols
    >>> m = symbols("m", integer=True, positive=True)
    >>> I_even = DirectSum(2*[L(m**2/2)], False)
    >>> I_even.dim == HC(m).dim
    True
    >>> [(I_even.rank - HC(m).rank).subs({m:k}) for k in [0,2,4,6,8,10]]
    [3, -3, 27, 237, 867, 2253]
    >>> I_odd  = DirectSum(2*[L((m**2 - 1)/2)] + [L(1)], False)
    >>> I_odd.dim == HC(m).dim
    True
    >>> [(I_odd.rank - HC(m).rank).subs({m:k}) for k in [1,3,5,7,9,11]]
    [2, -2, 86, 458, 1402, 3302]
    >>> I_three = L(9)
    >>> I_three.dim == HC(3).dim
    True
    >>> I_three.rank > HC(3).rank
    True

We can do the same for the real PSD cones, based on whether or not
``(m**2 + m) / 2`` is even or odd, again with a special case for ``m
== 3``::

    >>> I_even = DirectSum(2*[L((m**2 + m)/4)], False)
    >>> I_even.dim == fix_floor(HR(m).dim)
    True
    >>> [(I_even.rank - HR(m).rank).subs({m:k}) for k in [0,2,4,6,8,10]]
    [2, -2, 6, 64, 244, 630]
    >>> I_odd = DirectSum(2*[L( ((m**2 + m)/2 - 1)/2 )] + [L(1)], False)
    >>> I_odd.dim == fix_floor(HR(m).dim)
    True
    >>> [(I_odd.rank - HR(m).rank).subs({m:k}) for k in [1,3,5,7,9,11]]
    [2, -4, 20, 122, 384, 904]
    >>> I_three = L(6)
    >>> I_three.dim == HR(3).dim
    True
    >>> I_three.rank > HR(3).rank
    True

And the quaternion PSD cones::

    >>> I_even = DirectSum(2*[L((2*m**2 - m)/2)], False)
    >>> I_even.dim == fix_floor(HH(m).dim)
    True
    >>> [(I_even.rank - HH(m).rank).subs({m:k}) for k in [0,2,4,6,8,10]]
    [2, -8, 120, 914, 3286, 8532]
    >>> I_odd = DirectSum(2*[L((2*m**2 - m - 1)/2)] + [L(1)], False)
    >>> I_odd.dim == fix_floor(HH(m).dim)
    True
    >>> [(I_odd.rank - HH(m).rank).subs({m:k}) for k in [1,3,5,7,9,11]]
    [-1, 9, 365, 1787, 5379, 12629]

Finally, the 3x3 octonion cone::

   >>> I = DirectSum(3*[L(9)], False)
   >>> I.dim == HO(3).dim
   True
   >>> I.rank > HO(3).rank
   True


A systematic check of Theorem 5. For a given ``K.dim``,
``lowerbound3b`` is constant, ``lowerbound1`` is minimized by the
nonnegative orthant (with value ``n >= 2``), and ``lowerbound2`` is
minimized by the Lorentz cone with value ``K.dim + 2`` (proof: we are
either maximizing beta(K) or minimizing the Lyapunov rank in fixed
dimensions)::

   >>> from sympy import symbols, expand
   >>> d = symbols("d", integer=True, positive=True)
   >>> expand(fix_floor(lowerbound2(L(d))))
   d + 2

The second bound is increasing with ``K.dim`` and dominates the
first. As a result, we can use ``lowerbound2`` to determine the first
potentially valid ``n`` corresponding to any ``K.dim``. Moreover we
can easily compute the first ``K.dim`` where ``n + K.dim`` is
guaranteed to exceed the largest dimension we have cached by adding
``K.dim`` to the bound and setting it greater than that cached
dimension, i.e. by solving ``K.dim + 2 + K.dim > max_cone_dim()``::

    >>> from sql import all_cones_of_dim
    >>> result = True
    >>> max_d = (max_cone_dim() - 2) // 2
    >>> for d in range(1, max_d+1):
    ...     min_n = max(lowerbound1(RN(d)),
    ...                 lowerbound2(L(d)),
    ...                 lowerbound3b(L(d)))
    ...     max_n = max_cone_dim() - d
    ...     for n in range(min_n, max_n+1):
    ...         for K in all_cones_of_dim(d):
    ...             if n < lowerbound1(K): continue
    ...             if n < lowerbound2(K): continue
    ...             if n < lowerbound3b(K): continue
    ...             lhs = DirectSum([L(n),K])
    ...             for J in lhs.similacra():
    ...                 fs = list(J.factors())
    ...                 # remove() raises an error if L(n) isn't a factor!
    ...                 fs.remove(L(n))
    ...                 J_prime = DirectSum(fs)
    ...                 result &= J_prime in K.similacra()
    >>> result
    True

This is the function we'll use to check that most non-Lorentz factors
have a Lypaunov rank too small to be considered in Theorem 5. We need
to be careful that ``n + K.dim`` is at least large enough to hold the
cone, otherwise we might get false positives::

    >>> n,d = symbols("n,d", integer=True, positive=True)
    >>> f = lambda x: (x**2 - x + 2)/2   # no integer division... sympy
    >>> def rank_too_small(X):
    ...     g = f(n) + d - X.rank - f(n+d-X.dim)
    ...     return all( g.subs({n:i,d:j}) > 0
    ...                 for i in range(5,15)
    ...                 for j in range(1, max_dimK(i)+1)
    ...                 if i+j-X.dim >= 0 )

The argument to rule out ``HC(3)``, ``HR(4)``, ``HC(4)``, and
``HH(3)`` factors::

    >>> X = HC(3)
    >>> rank_too_small(X)
    True

    >>> X = HR(4)
    >>> rank_too_small(X)
    True

    >>> X = HC(4)
    >>> rank_too_small(X)
    True

    >>> X = HH(3)
    >>> rank_too_small(X)
    True

Also check the two cases that we argue via similacra::

    >>> DirectSum([HC(3),L(3),L(3)]) in HR(5).similacra()
    True
    >>> DirectSum([HC(3),L(4),L(4),L(3),L(1)]) in HR(6).similacra()
    True

Rule out multiple ``HR(3)`` factors to simplify the argument::

    >>> X = DirectSum(2*[HR(3)])
    >>> rank_too_small(X)
    True

And check that a single ``HR(3)`` is only viable as a factor of ``J``
when ``n >= 10`` and ``K.dim >= 6``::

    >>> X = HR(3)
    >>> g = f(n) + d - X.rank - f(n+d-X.dim)
    >>> all( g.subs({n:i,d:j}) > 0
    ...      for i in range(5,15)
    ...      for j in range(1, max_dimK(i)+1)
    ...      if (i < 10 or j < 6) )
    True

Check Theorem 5 directly by confirming that there simply aren't any
similacra for any of the valid cases::

    >>> from sql import all_cones_of_dim
    >>> all( not DirectSum([L(n),K]).similacra()
    ...      for n in range(5,15)
    ...      for d in range(1,max_dimK(n)+1)
    ...      for K in all_cones_of_dim(d)
    ...      if n >= max(lowerbound1(K),lowerbound2(K)) )
    True

Oh, and ensure we have enough cones computed for this to actually be
reliable::

    >>> from sql import max_cone_dim
    >>> 14 + max_dimK(14) <= max_cone_dim()
    True

Finally, we check the argument in Theorem 5 using the method that we
have described (partitions, possibly offset by one ``HR(3)``
factor). First, the small ``n`` where there are no ``HR(3)`` factors
to worry about. There are many matching signatures, but they're all
from isomorphic cones once you consider that ``L(2) == RN(2)``::

    >>> from signatures import partitions, f
    >>>
    >>> # We'll collect the matching signatures in a list
    >>> matches = []
    >>>
    >>> # reimplement the lower bounds in terms of partitions
    >>> lb1 = lambda p: 2 + partition_rank(p) - sum(p)
    >>> lb2 = lambda p: 2 + f(1+sum(p)) - partition_rank(p)
    >>> lb3 = lambda p: 5
    >>>
    >>> for n in range(5,10):
    ...     for d in range(1,max_dimK(n)+1):
    ...         for K in partitions(d):
    ...             if n < lb1(K): continue
    ...             if n < lb2(K): continue
    ...             if n < lb3(K): continue
    ...             Ln_K_rank = f(n) + partition_rank(K)
    ...             for J in partitions(d+n):
    ...                 if (partition_rank(J) == Ln_K_rank):
    ...                         Ln_K = [n] + K
    ...                         if not partitions_equivalent(J,Ln_K):
    ...                             matches.append( (Ln_K, J) )
    >>> matches
    []

Now things get a bit ugly, since we have to (potentially) include
``HR(3)`` in ``J``, which is always big enough to hold one. Whereas
before we represented a cone as a partition, we now represent it as an
``(j, p)`` pair, where ``j`` is either zero or one, indicating the
presence of an ``HR(3)`` factor, and ``p`` is a partition representing
its Lorentz factors::

    >>> matches = []
    >>> for n in range(10,15):
    ...     for d in range(1,max_dimK(n)+1):
    ...         Ks = []
    ...         for p in partitions(d):
    ...             # K is pure Lorentz, per the theorem
    ...             if n < lb1(p): continue
    ...             if n < lb2(p): continue
    ...             if n < lb3(p): continue
    ...             Ks.append( (0,p) )
    ...
    ...         # Always include pure-Lorentz J
    ...         Js = [ (0,p) for p in partitions(d+n) ]
    ...
    ...         # And since d+n is always >= 6, always include HR(3)
    ...         # with a partition of whatever's left.
    ...         for q in partitions(d + n - HR(3).dim):
    ...             Js.append( (1,q) )
    ...
    ...         for K in Ks:
    ...             Ln_K = (K[0], [n] + K[1])
    ...             Ln_K_rank = partition_rank(Ln_K[1])  # no HR(3)s here
    ...
    ...             for J in Js:
    ...                 J_rank = partition_rank(J[1]) + HR(3).rank*J[0]
    ...                 if J_rank == Ln_K_rank:
    ...                     # These two conditions aren't perfect, but
    ...                     # they're enough to eliminate all matches.
    ...                     if J[0] != Ln_K[0]:
    ...                         matches.append( (Ln_K,J) )
    ...                     if not partitions_equivalent(Ln_K[1],J[1]):
    ...                         matches.append( (Ln_K,J) )
    >>> matches
    []

"""

from signatures import *
from cones import *


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

    Being equivalent is a symmetric relationship::

        >>> from random import choice, randint
        >>> from signatures import partitions
        >>> n = randint(1,10)
        >>> ps = list(partitions(n))
        >>> p = choice(ps)
        >>> q = choice(ps)
        >>> partitions_equivalent(p,q) == partitions_equivalent(q,p)
        True

    """
    if (sum(p) != sum(q)):
        return False

    # sort and remove zeros (which shouldn't be there in the first
    # place)
    p = sorted(i for i in p if not i == 0)
    q = sorted(j for j in q if not j == 0)

    # Remove all factors of size >= 3 in p from both p and q
    for i in p:
        if i <= 2:
            continue
        elif i not in q:
            # Since i >= 3, if i is missing from q, they're not
            # isomorphic. (Without i >= 3 this doesn't work, because
            # for example [1,1] and [2] have no elements in common.)
            return False
        else:
            p.remove(i)
            q.remove(i)

    # Now what's left in p is its 1,2 elements; and what's left in q
    # is whatever 1,2 elements it had plus any elements >= 3 that
    # were not in p. We can repeat in the opposite direction.
    for j in q:
        if j <= 2:
            continue
        elif j not in p:
            return False
        else:
            p.remove(j)
            q.remove(j)

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


def max_dimK(n):
    r"""
    The largest possible dimension for ``K`` when ``n`` is fixed
    and must satisfy :func:`lowerbound1` and :func:`lowerbound2`.

    This is computed the dumb way rather than by solving the
    quadratic.

    Parameters
    ----------

    n : int
      The fixed dimension of the Lorentz cone in Theorem 5.

    Returns
    -------

    The largest integer dimension that ``K`` can have if ``n`` will
    satisfy both ``lowerbound1(K)`` and ``lowerbound2(K)``.

    Examples
    --------

    Check the table in Theorem 5::

        >>> max_dimK(5)
        3
        >>> max_dimK(6)
        4
        >>> max_dimK(7)
        4
        >>> max_dimK(8)
        5
        >>> max_dimK(9)
        5
        >>> max_dimK(10)
        6
        >>> max_dimK(11)
        6
        >>> max_dimK(12)
        6
        >>> max_dimK(13)
        7
        >>> max_dimK(14)
        7

    Compare the answer against the smart way, by solving the quadratic
    inequality ``d**2 - d - 4*n + 10 <= 0`` to find where the
    upwards-facing parabola last crosses the x-axis. (Note that
    we need ``n >= 3`` to get real solutions to this.)

        >>> from math import floor, sqrt
        >>> def qf(n):
        ...     a = 1
        ...     b = -1
        ...     c = 10 - 4*n
        ...     if (b**2 - 4*a*c) < 0:
        ...         return None
        ...     else:
        ...         return floor( (-b + sqrt(b**2 - 4*a*c)) / (2*a) )
        >>> all( qf(n) == max_dimK(n) for n in range(3,15) )
        True

    """
    d = 1
    while d*(d-1) <= (4*n - 10):
        d += 1
    return d-1


def partition_rank(p):
    r"""
    Return the Lyapunov rank of a sum of Lorentz factors whose
    dimensions are given by an integer partition.

    A direct sum of Lorentz cones is determined (almost) uniquely by
    the dimensions of its factors. The "almost" is because ``L(2)``
    and ``L(1) + L(1)`` are equal, but ``[1,1]`` and ``[2]`` are not.
    In any case, if we are given an integer partition that is intended
    to identify a direct sum of Lorentz cone (for example, computed by
    the :func:`signatures.partitions` function), then this function
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
    from signatures import f
    return sum( map(f,p) )


def fix_floor(s):
    r"""
    Strip symbolic "floor" calls from SymPy expressions.

    Sympy doesn't know (for example) that ``n**2 + n`` is even, so it
    will insert a :func:`sympy.functions.elementary.integers.floor`
    around the result of ``(n**2 + n)//2``. We however know that the
    "floor" is superfluous, so we can remove it using this function.

    Parameters
    ----------

    s : sympy.core.expr.Expr
      A sympy expression possibly containing a call to
      :func:`sympy.functions.elementary.integers.floor`.

    Returns
    -------

    Another :class:`sympy.core.expr.Expr`, devoid of floors.

    Examples
    --------

    ``n**2 + n`` is even::

        >>> from sympy import symbols
        >>> n = symbols("n", integer=True, positive=True)
        >>> expr = (n**2 + n) // 2
        >>> expr
        floor(n**2/2 + n/2)
        >>> fix_floor(expr)
        n**2/2 + n/2

    """
    from sympy.functions.elementary.integers import floor
    return s.replace(floor, lambda x: x)


def lowerbound1(K):
    r"""
    The first of the three lower bounds on "n" in Lemma 4, needed
    for Lemma 5 to hold.

    Parameters
    ----------

    K : SymmetricCone
      The cone for which you want the lower bound.

    Returns
    -------

    int
      The lower bound on "n" corresponding to ``K``.

    Examples
    --------

    Check Lemma 5 using our precomputed dict of admissible Lyapunov
    ranks. We start with a random cone ``K``, then add a Lorentz cone
    factor to it whose dimension is bounded below by this
    function. Then we loop through several values of "k", and check
    that none of them give rise to similacra: we subtract the
    signature of the L^{n+k} factor from both sides, and then check
    the precomputed list of signatures for the signature of what
    remains; basically, we exhaustively search for ``J`` in the Lemma.

        >>> from sql import (admissible_lorentz_ranks as alr,
        ...                  max_lorentz_rank_dim )
        >>> K = random_cone()
        >>> n_min = lowerbound1(K)
        >>> n_max = max_lorentz_rank_dim() - K.dim
        >>> k_max = min(K.dim, 20)
        >>> results = []
        >>> for n in range(n_min, n_max+1):
        ...     Lnplus = L(n)
        ...     lhs = DirectSum([Lnplus,K])
        ...     for k in range(1,k_max+1):
        ...         Lmplus = L(n+k)
        ...         target_dim = lhs.dim - Lmplus.dim
        ...         target_rank = lhs.rank - Lmplus.rank
        ...         results.append(target_rank not in alr(target_dim))
        >>> all(results)
        True

    This lower bound is tight, regardless of the other two::

        >>> K = random_cone()
        >>> while K.dim == 0:
        ...     K = random_cone()
        >>> n = lowerbound1(K) - 1
        >>> J1 = DirectSum([L(n), K])
        >>> J2 = DirectSum([L(n+1), RN(K.dim - 1)])
        >>> J1.signature() == J2.signature()
        True

    But we know some examples where this is the bound that's tight,
    and where ``K`` is big enough to ensure non-isomorphism (the
    existence of a real similacrum, not just a matching signature)::

        >>> m = 5
        >>> K = L(m)
        >>> n = lowerbound1(K) - 1
        >>> n >= lowerbound2(K)
        True
        >>> n >= lowerbound3b(K)
        True
        >>> J1 = DirectSum([L(n), K])
        >>> J2 = DirectSum([L(n+1), RN(K.dim - 1)])
        >>> J1.signature() == J2.signature()
        True
        >>> J1 == J2
        False

    We compute this lower bound in Proposition 11::

       >>> from sympy import symbols
       >>> m = symbols("m", integer=True, positive=True)
       >>> K = L(m)
       >>> fix_floor(lowerbound1(K)) == (m**2 - 3*m + 6)/2
       True

    """
    return 2 + K.rank - K.dim

def lowerbound2(K):
    r"""
    The second lower bound for "n" in Lemma 4.

    Parameters
    ----------

    K : SymmetricCone
      The cone for which you want the lower bound.

    Returns
    -------

    int
      The lower bound on "n" corresponding to ``K``.

    Examples
    --------

    It works with symbolic values, though we have to manually
    eliminate the `floor` that arises from Python's integer division
    (``m**2 + m + 2`` is guaranteed to be even)::

        >>> from sympy import expand, symbols
        >>> m = symbols("m", integer=True, positive=True)
        >>> expand(fix_floor(lowerbound2(L(m))))
        m + 2

    This bound is tight, because we have already computed
    a counterexample wherein the other two lower bounds are
    satisfied::

        >>> K = HC(3)
        >>> n = lowerbound2(K) - 1
        >>> n >= lowerbound1(K)
        True
        >>> n >= lowerbound3b(K)
        True
        >>> J1 = DirectSum([K,L(n)])
        >>> J2 = DirectSum([L(29),L(10)])
        >>> J1.signature() == J2.signature()
        True
        >>> J1 == J2
        False

    """
    from signatures import f
    return 2 + f(1 + K.dim) - K.rank


def lowerbound3a(K : SymmetricCone) -> int:
    r"""
    The third precondition on ``n`` in Lemma 4.

    Parameters
    ----------

    K : SymmetricCone
      The cone for which you want the lower bound (this parameter is
      essentially ignored).

    Returns
    -------

    int
      This lower bound is always 15.

    Examples
    --------

    Along with with the other two, this one implies a strict lower
    bound of ``2*K.dim``. If ``n`` satisfies the other two, we can add
    them up and divide by two to get another lower bound on ``n`` in
    terms of ``K.dim``. The assumption that ``n >= 15`` takes care of
    cones of dimension seven or less, and then the newly-derived
    inequality can be set greater than ``2*K.dim`` and solved to
    obtain a quadratic inequality that will be true for ``K.dim() >=
    8``.  Below we let the symbols ``d`` and ``r`` stand for ``K.dim``
    and ``K.rank``::

        >>> from sympy import expand, floor, symbols
        >>> d,r = symbols("d,r", integer=True, positive=True)
        >>> K = SymmetricCone(0)
        >>> K.dim = d
        >>> K.rank = r
        >>> implied_bound = 2*d
        >>> g = (lowerbound1(K) + lowerbound2(K))/2 - implied_bound

    We want ``g`` to be positive, but we can multiply it by ``4``
    without changing when it is positive. Again we have to strip
    the symbolic `floor` ourselves because sympy doesn't know that
    ``d**2 + d`` is even::

        >>> g = 4*g
        >>> fix_floor(expand(g))
        d**2 - 9*d + 10

    Since ``g`` is an upwards-facing parabola, it will be nonpositive
    on an interval, and positive everywhere else. We see that for ``d
    >= 8``, g will be positive::

        >>> [ g.subs({d:i}) for i in range(10) ]
        [10, 2, -4, -8, -10, -10, -8, -4, 2, 10]

    """
    return 15


def lowerbound3b(K):
    r"""
    The fourth and final precondition on ``n`` that we can use
    for Theorem 4, obtained near the end of the paper by loosening
    :func:`lowerbound3a`.

    Parameters
    ----------

    K : SymmetricCone
      The cone for which you want the lower bound (this parameter is
      essentially ignored).

    Returns
    -------

    int
      This lower bound is always 5.

    Examples
    --------

    Test the claim that if we drop the ``n >= 10`` condition, then the
    only counterexample we get is HR(3) ~ L(4) + L(2). We have to
    assume that ``n >= 3``, because that is implied by the first two
    bounds (just add them), but for ``n == 4`` we do get one
    counterexample. Thus the real bound is ``n >= 5``, for which we
    get no counterexamples::

        >>> from sql import all_cones_of_dim
        >>> from math import floor
        >>>
        >>> def check(n_start):
        ...     winners = []
        ...     for n in range(n_start, 10):
        ...         for d in range(1, floor(sqrt(4*n - 10))+2):
        ...             for K in all_cones_of_dim(d):
        ...                 if n < lowerbound1(K) or n < lowerbound2(K):
        ...                     continue
        ...                 C = DirectSum([L(n), K])
        ...                 for s in C.similacra():
        ...                     factors = [s]
        ...                     if isinstance(s,DirectSum):
        ...                         factors = s.factors()
        ...                     if L(n) not in factors:
        ...                         winners.append((K,n,s))
        ...     return winners
        >>>
        >>> check(3)
        [(L(1) + L(1), 4, HR(3))]
        >>> check(4)
        [(L(1) + L(1), 4, HR(3))]
        >>> check(5)
        []

    This bound cannot be lowered independent of the other two, which
    are both satisfied for the counterexample HR(3) ~ L(4) + L(2)::

        >>> K = RN(2)
        >>> n = lowerbound3b(K) - 1
        >>> n >= lowerbound1(K)
        True
        >>> n >= lowerbound2(K)
        True
        >>> DirectSum([K,L(n)]).similacra()
        (HR(3),)

    """
    return 5
