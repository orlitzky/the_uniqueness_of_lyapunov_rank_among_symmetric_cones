r"""
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
        >>> partitions_equivalent([1,3,3], [2,2,3])
        False

    Being equivalent is a symmetric relationship::

        >>> from random import randint
        >>> from signatures import partitions
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
    The first of the three lower bounds on "n", used in the
    proof of Lemma 3 and elsewhere.

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

    This lower bound is tight, regardless of the other two::

        >>> K = random_cone()
        >>> while K.dim == 0:
        ...     K = random_cone()
        >>> n = lowerbound1(K) - 1
        >>> J1 = DirectSum([L(n), K])
        >>> J2 = DirectSum([L(n+1), RN(K.dim - 1)])
        >>> J1.signature() == J2.signature()
        True

    We know some examples where this is the bound that's tight, and
    where ``K`` is big enough to ensure non-isomorphism (the existence
    of a real similacrum, not just a matching signature)::

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

    """
    return 2 + K.rank - K.dim


def lowerbound2(K):
    r"""
    The second lower bound on "n" used in the proof of Lemma 4
    and elsewhere.

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
    The third lower bound on "n", used in the proof of Lemma 2
    and elsewhere (and later loosened to :func:`lowerbound3b`).

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

    Yup::

        >>> lowerbound3a(random_cone())
        15

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


def test_proposition2():
    r"""
    Test the statement of Proposition 2.

    We construct a random cone and checking that its Lorentz rank is
    exceeded by that of the Lorentz cone of the same dimension.

    Setup
    -----

    Create some symbols that we'll use in the following examples::

        >>> from sympy import expand, simplify, symbols
        >>> m,n,k = symbols("m,n,k", integer=True, positive=True)

    Examples
    --------

    The result should hold::

        >>> test_proposition2()
        True

    Confirm that using two Lorentz cones gives you a smaller Lyapunov
    rank than one Lorentz cone would, assuming the dimension is
    fixed::

        >>> K = L(n)
        >>> J1 = L(k)
        >>> J2 = L(n-k)
        >>> J = DirectSum([J1,J2], False)  # can't sort symbolics
        >>> actual = fix_floor(K.rank - J.rank)
        >>> expected = (n-k)*k - 1
        >>> simplify(actual - expected)
        0

    The rank of ``L(27)`` exceeds that of ``HO(3)``::

        >>> K = L(27)
        >>> J = HO(3)
        >>> K.dim == J.dim
        True
        >>> K.rank > J.rank
        True

    The difference between the rank of the Lorentz cone and the rank
    of the real symmetric PSD cone of equal dimension is the
    polynomial we expect::

        >>> K = L((m**2 + m)/2)
        >>> J = HR(m)
        >>> fix_floor(K.dim - J.dim)
        0
        >>> expand(fix_floor(K.rank - J.rank))
        m**4/8 + m**3/4 - 9*m**2/8 - m/4 + 1

    The difference between the rank of the Lorentz cone and the rank
    of the complex Hermitian PSD cone of equal dimension is the
    polynomial we expect::

        >>> K = L(m**2)
        >>> J = HC(m)
        >>> fix_floor(K.dim - J.dim)
        0
        >>> fix_floor(K.rank - J.rank)
        m**4/2 - 5*m**2/2 + 2

    The difference between the rank of the Lorentz cone and the rank
    of the quaternion Hermitian PSD cone of equal dimension is the
    polynomial we expect::

        >>> K = L(2*m**2 - m)
        >>> J = HH(m)
        >>> fix_floor(K.dim - J.dim)
        0
        >>> expand(fix_floor(K.rank - J.rank))
        2*m**4 - 2*m**3 - 9*m**2/2 + m/2 + 1

    """
    K = random_cone()
    J = L(K.dim)

    # They could be isomorphic; but if not, it's a strict inequality.
    return ( J == K or J.rank > K.rank )


def test_proposition3() -> bool:
    r"""
    Test the statement of Proposition 3.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The result should hold::

        >>> test_proposition3()
        True

    From the proof of Proposition 3, we know that for ``n >= 3``,
    ``HR(n)`` has symmetric similacra with only Lorentz
    factors. Moreover when ``n < 3``, ``HR(n)`` _is_ a Lorentz
    cone. In either case, ``HR(n)`` should share its signature with a
    sum of Lorentz cones. We begin by computing the largest ``n`` for
    which we have the corresponding sum-of-Lorentz-cone data cached::

        >>> from sql import admissible_lorentz_ranks, max_lorentz_rank_dim
        >>> n_max = -1
        >>> if (mlrd := max_lorentz_rank_dim()) is not None:
        ...     while HR(n_max).dim <= mlrd:
        ...         n_max += 1
        ...     n_max -= 1
        >>> all(
        ...   HR(n).rank
        ...   in admissible_lorentz_ranks(HR(n).dim)
        ...   for n in range(n_max+1)
        ... )
        True

    In fact, we know the formula for at least one such similacrum::

        >>> def check(n):
        ...     K1 = L(n+1)
        ...     K2 = RN((n**2 - n - 2) // 2)
        ...     K = DirectSum([K1,K2])
        ...     return (HR(n).signature() == K.signature())
        >>>
        >>> all( check(n) for n in range(2,100) )
        True

    """
    from sql import max_cone_dim

    # Figure out how big "n" can be if we want to use the database of
    # cached cones (the ``similacra`` method uses it implicitly).
    n_max = 0
    if (mcd := max_cone_dim()) is not None:
        while HR(n_max).dim <= mcd:
            n_max += 1
        n_max -= 1

    n_min = 3
    if n_max <= n_min:
        # No cached cones?
        return True

    return all( HR(n).similacra() for n in range(n_min, n_max+1) )



def test_proposition4() -> bool:
    r"""
    Test the statement of Proposition 4.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The result should hold::

        >>> test_proposition4()
        True

    From the proof of Proposition 4, we know that for ``n >= 4``,
    ``HC(n)`` has symmetric similacra with only Lorentz
    factors. Moreover when ``n < 3``, ``HC(n)`` _is_ a Lorentz
    cone. In either case, ``HC(n)`` should share its signature with a
    sum of Lorentz cones. We begin by computing the largest ``n`` for
    which we have the corresponding sum-of-Lorentz-cone data cached
    (it's easy in this case)::

        >>> from math import floor, sqrt
        >>> from sql import admissible_lorentz_ranks, max_lorentz_rank_dim
        >>> n_max = -1
        >>> if (mlrd := max_lorentz_rank_dim()) is not None:
        ...     n_max = floor(sqrt(mlrd))
        >>> all(
        ...   HC(n).rank
        ...   in admissible_lorentz_ranks(HC(n).dim)
        ...   for n in range(n_max+1)
        ...   if n != 3
        ... )
        True

    We know the similacra for ``n >= 4`` explicitly; they are given in
    the proof of the proposition. There's a special case for ``n = 4``
    and then we handle ``n >= 5`` generically::

        >>> K = DirectSum(2*[L(5)] + [L(4), RN(2)])
        >>> K.signature() == HC(4).signature()
        True

        >>> all(
        ...   HC(n).signature() == K.signature()
        ...   for n in range(5,100)
        ...   if (K2 := RN(n**2 - 5*n + 1))
        ...   and (K := DirectSum(2*[L(n+1)] + [K2] + (n-1)*[L(3)]))
        ... )
        True

    An example using the :meth:`SymmetricCone.similacra` method for
    ``n = 4``::

        >>> from sql import max_cone_dim
        >>> K = DirectSum([L(5), L(5), L(4), RN(2)])
        >>> mcd = max_cone_dim()
        >>> skip = (not mcd) or (HC(4).dim > mcd)
        >>> skip or K in HC(4).similacra()
        True

    """
    from math import floor, sqrt
    from sql import max_cone_dim

    # Figure out how big "n" can be if we want to use the database of
    # cached cones (the ``similacra`` method uses it implicitly).
    n_max = 0
    if (mcd := max_cone_dim()) is not None:
        n_max = floor(sqrt(mcd))

    n_min = 4
    if n_max <= n_min:
        # No cached cones?
        return True

    return ( not HC(3).similacra()
             and
             all(HC(n).similacra() for n in range(n_min, n_max+1)) )



def test_proposition5() -> bool:
    r"""
    Test the statement of Proposition 5.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The result should hold::

        >>> test_proposition5()
        True

    From the proof of Proposition 5, we know that for ``n >= 3``,
    ``HH(n)`` has symmetric similacra with only Lorentz factors (the
    convenient ``HC(n+1)`` factor can be replaced by Lorentz cones
    using Proposition 4). Moreover when ``n < 3``, ``HH(n)`` _is_ a
    Lorentz cone. In either case, ``HH(n)`` should share its signature
    with a sum of Lorentz cones. We begin by computing the largest
    ``n`` for which we have the corresponding sum-of-Lorentz-cone data
    cached::

        >>> from sql import admissible_lorentz_ranks, max_lorentz_rank_dim
        >>> n_max = -1
        >>> if (mlrd := max_lorentz_rank_dim()) is not None:
        ...     while HH(n_max).dim <= mlrd:
        ...         n_max += 1
        ...     n_max -= 1
        >>>
        >>> all(
        ...   HH(n).rank
        ...   in admissible_lorentz_ranks(HH(n).dim)
        ...   for n in range(n_max+1)
        ... )
        True

    We know the formula for similacra explicitly; they are given in
    the proof of the proposition. There are special cases for ``n in
    [3,4,5]`` and then we handle ``n >= 6`` generically::

        >>> K = DirectSum([L(8), RN(7)])
        >>> K.signature() == HH(3).signature()
        True

        >>> K = DirectSum([L(10), RN(18)])
        >>> K.signature() == HH(4).signature()
        True

        >>> K = DirectSum([L(12), RN(33)])
        >>> K.signature() == HH(5).signature()
        True

        >>> all(
        ...   HH(n).signature() == K.signature()
        ...   for n in range(6,100)
        ...   if (K3 := RN(n**2 - 5*n - 3))
        ...   and (K := DirectSum([HC(n+1)] + 2*[L(n+1)] + [K3]))
        ... )
        True

    Further checks of the low-dimensional formulas using the
    :meth:`SymmetricCone.similacra` method::

        >>> from sql import max_cone_dim
        >>> mcd = max_cone_dim()

        >>> K = DirectSum([L(8), RN(7)])
        >>> skip = (not mcd) or (HH(3).dim > mcd)
        >>> skip or K in HH(3).similacra()
        True

        >>> K = DirectSum([L(10), RN(18)])
        >>> skip = (not mcd) or (HH(4).dim > mcd)
        >>> skip or K in HH(4).similacra()
        True

        >>> K = DirectSum([L(12), RN(33)])
        >>> skip = (not mcd) or (HH(5).dim > mcd)
        >>> skip or K in HH(5).similacra()
        True

    """
    from sql import max_cone_dim

    # Figure out how big "n" can be if we want to use the database of
    # cached cones (the ``similacra`` method uses it implicitly).
    n_max = 0
    if (mcd := max_cone_dim()) is not None:
        while HH(n_max).dim <= mcd:
            n_max += 1
        n_max -= 1

    n_min = 3
    if n_max <= n_min:
        # No cached cones?
        return True

    return all( HH(n).similacra() for n in range(n_min, n_max+1) )


def test_proposition6() -> bool:
    r"""
    Test the statement of Proposition 6.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The result should hold::

        >>> test_proposition6()
        True

    From the proof of Proposition 6, we know that ``HO(3)`` has a
    symmetric similacrum with only Lorentz factors. Moreover when ``n
    < 3``, ``HO(n)`` _is_ a Lorentz cone. In either case, ``HO(n)``
    should share its signature with a sum of Lorentz cones. We begin
    by computing the largest ``n`` for which we have the corresponding
    sum-of-Lorentz-cone data cached::

        >>> from sql import admissible_lorentz_ranks, max_lorentz_rank_dim
        >>> n_max = -1
        >>> if (mlrd := max_lorentz_rank_dim()) is not None:
        ...     while n_max <= 3 and HO(n_max).dim <= mlrd:
        ...         n_max += 1
        ...     n_max -= 1
        >>>
        >>> all(
        ...   HO(n).rank
        ...   in admissible_lorentz_ranks(HO(n).dim)
        ...   for n in range(n_max+1)
        ... )
        True

    We know the formula for one similacrum explicitly; it is given
    in the proof of the proposition::

        >>> K = DirectSum([L(11),L(5),L(3),RN(8)])
        >>> K.signature() == HO(3).signature()
        True

    Repeat with the cached similacra data::

        >>> from sql import have_cone_dim
        >>> (not have_cone_dim(HO(3).dim)) or K in HO(3).similacra()
        True

    """
    from sql import max_cone_dim

    if HO(3).dim > max_cone_dim():
        # no data
        return True

    return not (not HO(3).similacra())


def test_theorem2() -> bool:
    r"""
    Test the statement of Theorem 2.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The result should hold::

        >>> test_theorem2()
        True

    Implicit in the proof of this theorem is the fact that every
    irreducible symmetric cone other than ``HC(3)`` shares its
    signature with a sum or Lorentz cones::

        >>> from sql import admissible_lorentz_ranks, max_lorentz_rank_dim
        >>> Ks = ( random_irreducible_cone() for _ in range(100) )
        >>> all ( K.rank in admissible_lorentz_ranks(K.dim)
        ...       or K == HC(3)
        ...       for K in Ks
        ...       if K.dim <= max_lorentz_rank_dim() )
        True

    """
    from sql import max_cone_dim
    Ks = ( random_irreducible_cone() for _ in range(100) )
    return all ( not (not K.similacra())
                 or K == HC(3)
                 or isinstance(K,L)
                 for K in Ks
                 if K.dim <= max_cone_dim() )


def test_lemma1() -> bool:
    r"""
    Test the statement of Lemma 1.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The result should hold::

        >>> test_lemma1()
        True

    """
    from sql import max_cone_dim
    Ks = ( random_cone() for _ in range(100) )

    return all(
      K.similacra()
      or
      all(not K_i.similacra() for K_i in K.factors())
      for K in Ks
      if K.dim <= max_cone_dim()
    )


def test_corollary2() -> bool:
    r"""
    Test the statement of Corollary 2.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The result should hold::

        >>> test_corollary2()
        True

    """
    from sql import max_cone_dim

    Ks = ( random_cone() for _ in range(100) )
    return all(
      not (not K.similacra())
      for K in Ks
      if K.dim <= max_cone_dim()
      and K.factors().count(HC(3)) > 1
    )


def test_theorem3() -> bool:
    r"""
    Test the statement of Theorem 3.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The result should hold::

        >>> test_theorem3()
        True

    """
    from signatures import partitions

    # The fact that every cone shares a signature with a sum of
    # Lorentz cones and/or HC(3) is the basis for the function
    # sql.admissible_ranks(), but here we test it directly using
    # partitions.
    K = random_cone()
    while K.dim > 60:
        # Make sure we don't have to partition anything too big (it
        # takes a looong time). Partitioning 60 can be done in a few
        # seconds, but e.g. 80 may crash the machine.
        K = random_cone()

    r = False
    r |= any( partition_rank(p) == K.rank for p in partitions(K.dim) )

    dim_rest = K.dim - 9
    if dim_rest >= 0:
        # Try with an HC(3) factor
        r |= any( partition_rank(p) == (K.rank - 17)
                  for p in partitions(dim_rest) )
    return r


def test_example1() -> bool:
    r"""
    Test Example 1.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The relationships in the example should hold::

        >>> test_example1()
        True

    """
    result = True

    # The two similacra of HC(3) mentioned in the example
    K1 = DirectSum([L(11),L(3)] + [L(5),RN(8)])
    K2 = DirectSum([L(11),L(3)] + [L(4)] + 3*[L(3)])

    # Signature comparison, suffices because they're obviously
    # not isomorphic
    result &= ( K1.signature() == K2.signature() )
    result &= ( K1.signature() == HO(3).signature() )

    # And repeat using cached similacra if possible
    from sql import have_cone_dim
    if have_cone_dim(HO(3).dim):
        result &= K1 in K2.similacra()
        result &= K2 in K1.similacra()

        result &= K1 in HO(3).similacra()
        result &= HO(3) in K1.similacra()

        result &= K2 in HO(3).similacra()
        result &= HO(3) in K2.similacra()

    return result


def test_lemma2() -> bool:
    r"""
    Test Lemma 2.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The implication in the result should hold::

        >>> test_lemma2()
        True

    In the proof of this lemma, we average :func:`lowerbound1` and
    :func:`lowerbound2` to obtain a new lower bound on ``n``, and then
    set that new lower bound greater than ``2*K.dim``. This leads to a
    quadratic inequality (which we call ``g`` below) that can easily
    be solved, and will hold for ``K.dim >= 8``. The ``K.dim < 7``
    cases follow trivially from :func:`lowerbound3a`. Below we let the
    symbols ``d`` and ``r`` stand for ``K.dim`` and ``K.rank``::

        >>> from sympy import expand, floor, symbols
        >>> d,r = symbols("d,r", integer=True, positive=True)
        >>> K = SymmetricCone(0)
        >>> K.dim = d
        >>> K.rank = r
        >>> g = (lowerbound1(K) + lowerbound2(K))/2 - 2*d

    We want ``g`` to be positive, i.e. for the new average bound to be
    strictly greater than ``2*d``. We can multiply it by ``4`` without
    changing when it is positive. Again we have to strip the symbolic
    ``floor`` ourselves because sympy doesn't know that ``d**2 + d``
    is even::

        >>> g = 4*g
        >>> fix_floor(expand(g))
        d**2 - 9*d + 10

    Since ``g`` is an upwards-facing parabola, it will be nonpositive
    on an interval, and positive everywhere else. We see that for ``d
    >= 8``, ``g`` will be positive::

        >>> def sgn(x):
        ...     if x < 0: return -1
        ...     elif x == 0: return  0
        ...     else: return  1
        >>> for i in range(12):
        ...     print(f"d = {i : >2}, sgn(g(d)) = {sgn(g.subs({d:i})) : >2}")
        d =  0, sgn(g(d)) =  1
        d =  1, sgn(g(d)) =  1
        d =  2, sgn(g(d)) = -1
        d =  3, sgn(g(d)) = -1
        d =  4, sgn(g(d)) = -1
        d =  5, sgn(g(d)) = -1
        d =  6, sgn(g(d)) = -1
        d =  7, sgn(g(d)) = -1
        d =  8, sgn(g(d)) =  1
        d =  9, sgn(g(d)) =  1
        d = 10, sgn(g(d)) =  1
        d = 11, sgn(g(d)) =  1

    """
    from random import randint
    K_n_pairs = ( (random_cone(),randint(0,1000))
                  for _ in range(100) )
    return all(
      (n > 2*K.dim) or any([n < lowerbound1(K),
                            n < lowerbound2(K),
                            n < lowerbound3a(K)])
      for (K,n) in K_n_pairs
    )


def test_lemma3() -> bool:
    r"""
    Test Lemma 3.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The implication in the result should hold::

        >>> test_lemma3()
        True

    """
    from random import randint
    from sql import admissible_ranks, max_lorentz_rank_dim

    # Repeat the check 100 times with random values. If any of them
    # fail, we set the result to False before returning it.
    result = True
    Ks = ( random_cone() for _ in range(100) )

    for K in Ks:
        # The check for a single cone. We start with a random cone
        # ``K``, then add a Lorentz cone factor to it whose dimension
        # is bounded below according to the Lemma.
        n_min = lowerbound1(K)

        # We're going to use sql.admissible_ranks() to find signature
        # matches for K + L(n), so we have to make sure that we have
        # rank data cached up to dimension ``K.dim + n`` at least.
        n_max = max_lorentz_rank_dim() - K.dim
        if n_min > n_max:
            continue

        # There's room for an ``L(n)``, so let's add it.
        n = randint(n_min, n_max)
        Ln_plus_K = DirectSum([L(n), K])

        # The largest possible value of k is K.dim: if k exceeds
        # K.dim, then we wind up searching for a J whose dimension is
        # negative. The Lemma requires k >= 1, though, so we may skip
        # trivial K.
        if K.dim == 0:
            continue

        # Otherwise, pick a k at random.
        k = randint(1, K.dim)

        # The dimension and rank that J must have in the Lemma.
        J_dim = Ln_plus_K.dim - L(n+k).dim
        J_rank = Ln_plus_K.rank - L(n+k).rank

        # Does such a J exist? It should not.
        result &= J_rank not in admissible_ranks(J_dim)

    return result



def test_lemma4() -> bool:
    r"""
    Test Lemma 4.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The implication in the result should hold::

        >>> test_lemma4()
        True

    Check the symbolic identity for the derivative of the g-delta
    function in this Lemma::

        >>> from sympy import diff, symbols
        >>> x,d = symbols("x,d", integer=True, positive=True)
        >>> g = fix_floor(L(x).rank - L(x-d).rank)
        >>> diff(g, x)
        d

    Check inequality (1) with a few concrete examples::

        >>> from random import randint
        >>> from signatures import f
        >>> x = randint(0,30)
        >>> y = randint(0,30)
        >>> f(x+y) >= f(x) + f(y)
        True

    Check the symbolic identity (2) for the Lyapunov rank of
    ``L(n-1)`` in terms of that of ``L(n)`` in this Lemma::

        >>> from sympy import simplify, symbols
        >>> n = symbols("n", integer=True, positive=True)
        >>> lhs = L(n-1).rank
        >>> rhs = L(n).rank - (n-1)
        >>> simplify(fix_floor(lhs - rhs))
        0

    Check implication (3)::

        >>> f = lambda x: L(x).rank
        >>> all( f(n-1) + f(1+dimK) >= f(n-1-d) + f(1+dimK+d)
        ...      for d in range(100)
        ...      for dimK in range(100)
        ...      for n in range(2 + dimK + d, 100) )
        True

    """
    from random import randint
    from signatures import partitions

    Ks = ( random_cone() for _ in range(10) )

    result = True
    for K in Ks:
        # Partitioning 60 can be done in a few seconds, but e.g. 80
        # may crash the machine. Keep in mind that we're partitioning
        # n + dim(K), and that this is multiplied (at worst) by the
        # length of Ks!
        min_n = max(lowerbound1(K), lowerbound2(K), lowerbound3a(K))
        max_n = 60 - K.dim

        if min_n > max_n:
            # This is a pretty tight window. For example we know that
            # min_n = 31 for K == HC(3), but then K.dim == 9 already.
            continue
        n = randint(min_n, max_n)

        J = DirectSum([L(n),K])
        for p in partitions(n + K.dim, n-1):
            result &= (partition_rank(p) < J.rank)

    return result



def test_theorem4() -> bool:
    r"""
    Do nothing and return ``True``.

    Theorem 5 is a stronger version of Theorem 4, so there is no point
    in testing Theorem 4, on its own, directly. There are however some
    doctests for the proof strategy below.

    Examples
    --------

        >>> test_theorem4()
        True

    In the proof of this theorem, we "replace" the non-Lorentz
    irreducible factors with sums of Lorentz cones. Here we confirm
    that those sums have the correct dimensions, and Lyapunov ranks
    that dominates the Lyapunov ranks of the things they replace. The
    first example we give is for the complex PSD cones::

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
    ``(m**2 + m) / 2`` is even or odd, again with a special case for
    ``m == 3``::

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

    """
    return True


def test_theorem5() -> bool:
    r"""
    Test Theorem 5.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The implication in the result should hold::

        >>> test_theorem5()
        True

    For a given ``K.dim``, :func:`lowerbound3b` is constant,
    :func:`lowerbound1` is minimized by the nonnegative orthant (with
    value ``n >= 2``), and :func:`lowerbound2` is minimized by the
    Lorentz cone with value ``K.dim + 2`` (proof: we are either
    maximizing beta(K) or minimizing the Lyapunov rank in fixed
    dimensions). The second bound is therefore increasing with
    ``K.dim`` and obviously dominates the first. As a result, we can
    use :func:`lowerbound2` to determine the first potentially valid
    ``n`` corresponding to any ``K.dim``. Moreover we can use this to
    compute the first ``K.dim`` such that ``n + K.dim`` will exceed
    the largest dimension we have cached: basically we just add
    ``K.dim`` to the bound and set it greater than the largest
    dimension we have cached, i.e. we solve ``K.dim + 2 + K.dim >
    max_cone_dim()``::

        >>> lowerbound1(RN(3))
        2
        >>> lowerbound1(RN(8))
        2
        >>> lowerbound1(RN(22))
        2
        >>> lowerbound1(RN(57))
        2
        >>> from sympy import symbols, expand
        >>> d = symbols("d", integer=True, positive=True)
        >>> expand(fix_floor(lowerbound2(L(d))))
        d + 2

    Now we begin to verify the techniques used in the proof. First we
    define the function that we'll use to confirm that most
    non-Lorentz factors have a Lypaunov rank too small to be
    considered in the theorem. We need to be careful that ``n +
    K.dim`` is at least large enough to hold a particular factor,
    otherwise we might get false positives::

        >>> from sympy import symbols
        >>> n,d = symbols("n,d", integer=True, positive=True)
        >>> f = lambda x: (x**2 - x + 2)/2   # no integer division... sympy
        >>> def rank_too_small(X):
        ...     g = f(n) + d - X.rank - f(n+d-X.dim)
        ...     return all( g.subs({n:i,d:j}) > 0
        ...                 for i in range(5,15)
        ...                 for j in range(1, max_dimK(i)+1)
        ...                 if i+j-X.dim >= 0 )

    The argument we use to rule out ``HC(3)``, ``HR(4)``, ``HC(4)``,
    and ``HH(3)`` factors::

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

    The two cases that we argue via similacra::

        >>> from sql import max_cone_dim
        >>> ( HR(5).dim > max_cone_dim()
        ...   or
        ...   DirectSum([HC(3),L(3),L(3)]) in HR(5).similacra() )
        True
        >>> ( HR(6).dim > max_cone_dim()
        ...   or
        ...   DirectSum([HC(3),L(4),L(4),L(3),L(1)]) in HR(6).similacra() )
        True

    Rule out multiple ``HR(3)`` factors to simplify the argument::

        >>> X = DirectSum(2*[HR(3)])
        >>> rank_too_small(X)
        True

    A single ``HR(3)`` is only viable as a factor of ``J`` when ``n >=
    10`` and ``K.dim >= 6`` (to check this, we partially reimplement
    ``rank_too_small`` to restrict ``n`` and ``K.dim`` accordingly)::

        >>> X = HR(3)
        >>> g = f(n) + d - X.rank - f(n+d-X.dim)
        >>> all( g.subs({n:i,d:j}) > 0
        ...      for i in range(5,15)
        ...      for j in range(1, max_dimK(i)+1)
        ...      if (i < 10 or j < 6) )
        True

    One of the last statements in the proof is that the conclusion is
    easy to verify for the "new" cases because there simply aren't any
    new similacra (so we don't even have to worry about whether or not
    "J" has the stated form). We're checking that a list of similacra
    is empty, so there is no need to check the dimension of our cone
    against :func:`sql.max_cone_dim`; if we exceed it, we'll get back
    empty lists of similacra anyway::

        >>> from sql import all_cones_of_dim
        >>> all( not DirectSum([L(n),K]).similacra()
        ...      for n in range(5,15)
        ...      for d in range(1,max_dimK(n)+1)
        ...      for K in all_cones_of_dim(d)
        ...      if  n >= max(lowerbound1(K),lowerbound2(K)) )
        True

    Finally, we check the proof using the low-tech method that we have
    described: partitions, possibly offset by one ``HR(3)``
    factor. First, the small ``n`` where there are no ``HR(3)``
    factors to worry about. There are many matching signatures, but
    they're all from isomorphic cones once you consider that ``L(2) ==
    RN(2)``::

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
    ``HR(3)`` in ``J``, which is always big enough to hold
    one. Whereas before we represented a cone as a partition, we now
    represent it as an ``(j, p)`` pair, where ``j`` is either zero or
    one, indicating the presence of an ``HR(3)`` factor, and ``p`` is
    a partition representing its Lorentz factors::

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
        ...             Ln_K_rank = partition_rank(Ln_K[1])  # no HR(3) here
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

    Check the relationships between the lower bounds on ``n``. Each
    can be violated while the others are satisfied, and in each case,
    Theorem 5 fails. This is discussed subsequent to Theorem 5 in the
    paper::

        >>> from random import randint
        >>> m = randint(5,20)
        >>> K = L(m)
        >>> lowerbound3b(K) < lowerbound2(K) < lowerbound1(K)
        True
        >>> n = lowerbound1(K) - 1
        >>> (n >= lowerbound1(K), n >= lowerbound2(K), n >= lowerbound3b(K))
        (False, True, True)
        >>> ( DirectSum([L(n),L(m)]).signature()
        ...   ==
        ...   DirectSum([L(n+1),RN(m-1)]).signature() )
        True

        >>> K = RN(2)
        >>> n = 4
        >>> lowerbound1(K) < lowerbound2(K) < lowerbound3b(K)
        True
        >>> (n >= lowerbound1(K), n >= lowerbound2(K), n >= lowerbound3b(K))
        (True, True, False)
        >>> DirectSum([L(n),K]).signature() == HR(3).signature()
        True

        >>> K = HC(3)
        >>> lowerbound3b(K) < lowerbound1(K) < lowerbound2(K)
        True
        >>> n = 30
        >>> (n >= lowerbound1(K), n >= lowerbound2(K), n >= lowerbound3b(K))
        (True, False, True)
        >>> ( DirectSum([L(n),K]).signature()
        ...   ==
        ...   DirectSum([L(29),L(10)]).signature() )
        True

    """
    # Use cached data so that we can go beyond the n=15 case
    # and check the result for Theorem 4, too.
    from sql import all_cones_of_dim, max_cone_dim
    result = True

    # Explained in the docstring. We don't _really_ have to stop where
    # L(n)+K will have the max cached dim, because anything larger
    # will appear to have no similacra, and (for the sake of the
    # theorem) that's fine. But we do have to stop _somewhere_, so we
    # might as well stop here?
    max_d = (max_cone_dim() - 2) // 2

    for d in range(1, max_d+1):
        min_n = max(lowerbound1(RN(d)),
                    lowerbound2(L(d)),
                    lowerbound3b(L(d)))
        max_n = max_cone_dim() - d
        for n in range(min_n, max_n+1):
            for K in all_cones_of_dim(d):
                if n < lowerbound1(K): continue
                if n < lowerbound2(K): continue
                if n < lowerbound3b(K): continue
                lhs = DirectSum([L(n),K])
                for J in lhs.similacra():
                    fs = list(J.factors())
                    # remove() raises an error if L(n) isn't a factor,
                    # so this guarantees that L(n) is one.
                    fs.remove(L(n))
                    J_prime = DirectSum(fs)
                    result &= J_prime in K.similacra()

    return result


def test_corollary3() -> bool:
    r"""
    Test Corollary 3.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The implication in the result should hold::

        >>> test_corollary3()
        True

    Check the claim for ``HC(3)`` directly, using whatever cached
    cones are available::

        >>> from sql import max_cone_dim
        >>> n_min = 31
        >>> n_max = max_cone_dim() - 9
        >>> all( not DirectSum([L(n),HC(3)]).similacra()
        ...      for n in range(n_min, n_max+1) )
        True

    If you think hard about it, or consult an earlier version of the
    paper, you will conclude that only sums of Lorentz cones need to
    be checked for similacra of ``HC(3) + L(n)``. Here we repeat the
    check above using our cached Lorentz ranks (which are easier to
    compute)::

        >>> from sql import admissible_lorentz_ranks, max_lorentz_rank_dim
        >>> n_min = 31
        >>> n_max = max_lorentz_rank_dim() - 9
        >>>
        >>> all(
        ...   (17+L(n).rank)
        ...   not in admissible_lorentz_ranks(9+n)
        ...   for n in range(n_min, n_max+1)
        ... )
        True

    """
    from random import randint
    from sql import all_cones_of_dim, max_cone_dim

    # We'll logical-and this with the result from each test case.
    result = True

    # Use the same trick we used in test_theorem5() to decide where to
    # stop.
    max_d = (max_cone_dim() - 2) // 2

    for d in range(1, max_d+1):
        max_n = max_cone_dim() - d  # leave room for K
        for K in all_cones_of_dim(d):
            min_n = max(lowerbound1(K),lowerbound2(K),lowerbound3b(K))
            if min_n > max_n:
                # The dimension of L(n)+K is guaranteed to exceed the
                # dimensions we have cached, so skip this cone.
                continue

            if not K.similacra():
                # the corollary says that L(n)+K should have no
                # similacra
                n = randint(min_n, max_n)
                result &= not DirectSum([L(n), K]).similacra()

    return result


def test_proposition7() -> bool:
    r"""
    Test Proposition 7.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The implication in the result should hold::

        >>> test_proposition7()
        True

    Confirm the table for ``n <= 30`` using sums of Lorentz cones. If
    you think hard enough about it (or read an earlier version of the
    paper), this suffices::

        >>> from sql import admissible_lorentz_ranks, max_lorentz_rank_dim
        >>>
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
        ...         K = DirectSum([ L(i) for i in d[n] ])
        ...         return(
        ...           J.rank in admissible_lorentz_ranks(J.dim)
        ...           and
        ...           K.signature() == J.signature()
        ...         )
        ...     else:
        ...         return J.rank not in admissible_lorentz_ranks(J.dim)
        >>>
        >>> n_max = max_lorentz_rank_dim() - 9
        >>> all( check(n) for n in range(n_max+1) )
        True

    """
    from sql import max_cone_dim

    # the list of "n" where we expect to find similacra
    expected_n = [2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 18, 21, 22, 30]

    # We're going to be calling similacra() on L(n)+HC(3), so make
    # sure we don't make "n" so large that we exceed what is cached.
    n_max = max_cone_dim() - 9

    return all(
      (not DirectSum([HC(3),L(n)]).similacra())
      or
      (n in expected_n)
      for n in range(n_max+1)
    )


def test_proposition8() -> bool:
    r"""
    Test Proposition 8.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The implication in the result should hold::

        >>> test_proposition8()
        True

    Verify the cases mentioned explicitly in the proof. First, the
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

    And now the ``m == 2`` case where there is exactly one
    counterexample::

        >>> m = 2
        >>> K = L(m)
        >>> lowerbound1(K)
        2
        >>> DirectSum([K,L(2)]).similacra()
        ()
        >>> DirectSum([K,L(3)]).similacra()
        ()
        >>> from sql import max_cone_dim
        >>> expected = (HR(3),)
        >>> J = DirectSum([K,L(4)])
        >>> (J.dim > max_cone_dim()) or (J.similacra() == expected)
        True

   Confirm that the stated bound actually comes from
   :func:`lowerbound1`::

       >>> from sympy import symbols
       >>> m = symbols("m", integer=True, positive=True)
       >>> fix_floor(lowerbound1(L(m))) == (m**2 - 3*m + 6)/2
       True

    """
    from sql import max_cone_dim
    mcd = max_cone_dim()

    return all(
      not DirectSum([L(m),L(n)]).similacra()
      for m in range(1, mcd + 1)
      for n in range(lowerbound1(L(m)), mcd - m)
      if m != 2
    )
