from partitions import *
from cones import *


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

        >>> max_dimK(3)
        2
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
    proof of Lemma 4 and elsewhere.

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
        >>> J1 = DirectSum([K, L(n)])
        >>> J2 = DirectSum([RN(K.dim - 1), L(n+1)])
        >>> J1.signature() == J2.signature()
        True

    We know some examples where this is the bound that's tight, and
    where ``K`` is big enough to ensure non-isomorphism (the existence
    of a real simulacrum, not just a matching signature)::

        >>> m = 5
        >>> K = L(m)
        >>> n = lowerbound1(K) - 1
        >>> n >= lowerbound2(K)
        True
        >>> n != 4
        True
        >>> J1 = DirectSum([K, L(n)])
        >>> J2 = DirectSum([RN(K.dim - 1), L(n+1)])
        >>> J1.signature() == J2.signature()
        True
        >>> J1 == J2
        False

    """
    return 2 + K.rank - K.dim


def lowerbound2(K):
    r"""
    The second lower bound on "n" used in the proof of Lemma 5
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
        >>> n != 4
        True
        >>> J1 = DirectSum([K,L(n)])
        >>> J2 = DirectSum([L(29),L(10)])
        >>> J1.signature() == J2.signature()
        True
        >>> J1 == J2
        False

    This lower bound is exactly what is needed in Lemma 5::

        >>> from sympy import expand, symbols
        >>> n,d,r = symbols("n,d,r", integer=True, positive=True)
        >>> K = SymmetricCone(0)
        >>> K.dim = d
        >>> K.rank = r
        >>> ineq1 = f(n) - f(n-1) + 1 >= 2 + f(K.dim + 1) - K.rank
        >>> ineq1 = expand(fix_floor(ineq1))
        >>> ineq1 == expand(fix_floor(n >= lowerbound2(K)))
        True

    """
    return 2 + f(1 + K.dim) - K.rank


def lowerbound3(K : SymmetricCone) -> int:
    r"""
    The third lower bound on "n", used in the proof of Lemma 3
    and elsewhere (and later loosened to ``n != 4``).

    Parameters
    ----------

    K : SymmetricCone
      The cone for which you want the lower bound (this parameter is
      essentially ignored).

    Returns
    -------

    int
      This lower bound is always 10.

    Examples
    --------

    Yup::

        >>> lowerbound3(random_cone())
        10

    """
    return 10


def test_lemma1():
    r"""
    Test the statement of Lemma 1.

    We construct many random irreducible cones and checking that their
    signatures match if and only if they are equal.

    Examples
    --------

    The result should hold::

        >>> test_lemma1()
        True

    The only cone of the same dimension as ``HO(3)`` is ``L(27)``, and
    its Lyapunov rank is too large::

        >>> HR(6).dim, HR(7).dim
        (21, 28)
        >>> HC(5).dim, HC(6).dim
        (25, 36)
        >>> HH(3).dim, HH(4).dim
        (15, 28)
        >>> L(27).rank
        352

    There are no members (``3``-by-``3`` or larger) of the matrix
    families with matching signatures::

        >>> from sympy import symbols, solve
        >>> m,n = symbols("m,n", integer=True, positive=True)
        >>> solve([fix_floor(HR(m).dim) - HC(n).dim, HR(m).rank - HC(n).rank], n)
        []
        >>> solve([fix_floor(HR(m).dim) - HC(n).dim, HR(m).rank - HC(n).rank], m)
        []
        >>> solve([fix_floor(HR(m).dim) - HH(n).dim, HR(m).rank - HH(n).rank], n)
        []
        >>> solve([fix_floor(HR(m).dim) - HH(n).dim, HR(m).rank - HH(n).rank], m)
        []
        >>> solve([HC(m).dim - HH(n).dim, HC(m).rank - HH(n).rank], m)
        []
        >>> solve([HC(m).dim - HH(n).dim, HC(m).rank - HH(n).rank], n)
        []

    We should actually be obtaining the solution ``m == n == 1`` for
    all of these, but the formula (per Appendix A) for the Lyapunov
    rank of ``HH(1)`` is wrong.  I don't know why the solution set is
    empty for ``HR(1)`` and ``HC(1)``. We can easily verify::

        >>> [HR(1).dim - HC(1).dim, HR(1).rank - HC(1).rank]
        [0, 0]

    If we allow ``m`` and ``n`` to be zero, the ``m == n == 1``
    solution shows up indirectly::

        >>> m,n = symbols("m,n", integer=True, nonnegative=True)
        >>> solve([HR(m).dim - HC(n).dim, HR(m).rank - HC(n).rank], m)
        [(-sqrt(2*n**2 - 1),), (sqrt(2*n**2 - 1),)]

    However, only the latter of these is real and positive, and only
    for ``n >= 1``. If we substitute ``m == sqrt(2*n**2 - 1)`` into
    the equation relating the dimensions, we derive ``n == 1``, albeit
    with no help from SymPy, who thinks that the equation
    ``sqrt(2*n**2 - 1) == 1`` has no solutions::

        >>> from sympy import sqrt
        >>> eq = 2*(fix_floor(HR(sqrt(2*n**2 - 1)).dim) - HC(n).dim)
        >>> eq
        sqrt(2*n**2 - 1) - 1
        >>> solve(eq, n)
        []

    Anyway, back to the argument: the signature of a Lorentz cone
    matches that of an ``n``-by-``n`` matrix cone only for ``n <= 2``,
    where they are isomorphic to Lorentz cones. Though note that the
    formulas for Lyapunov rank (per Appendix A) are not even valid for
    the matrix cones when ``n <= 1``::

        >>> solve(fix_floor(L(HR(n).dim).rank - HR(n).rank), n)
        [1, 2]
        >>> solve(fix_floor(L(HC(n).dim).rank - HC(n).rank), n)
        [1, 2]
        >>> solve(fix_floor(L(HH(n).dim).rank - HH(n).rank), n)
        [2]

    """
    # Use lists here to avoid quietly exhausting the generator in the
    # "all" comprehension.
    Ks = [ random_irreducible_cone(0,20) for _ in range(250) ]
    Js = [ random_irreducible_cone(0,20) for _ in range(250) ]

    return all(
      (K.signature() == J.signature())
      ==
      (K == J)
      for K in Ks
      for J in Js
    )


def test_proposition3():
    r"""
    Test the statement of Proposition 3.

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

        >>> test_proposition3()
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

    From the proof of Proposition 4, we know that for ``n >= 3``,
    ``HR(n)`` has symmetric simulacra with only Lorentz
    factors. Moreover when ``n < 3``, ``HR(n)`` _is_ a Lorentz
    cone. In either case, ``HR(n)`` should share its signature with a
    sum of Lorentz cones. We begin by computing the largest ``n`` for
    which we have the corresponding sum-of-Lorentz-cone data cached::

        >>> from sql import admissible_lorentz_ranks, max_lorentz_rank_dim
        >>> n_max = -1
        >>> mlrd = max_lorentz_rank_dim()
        >>> while HR(n_max).dim <= mlrd:
        ...     n_max += 1
        >>> n_max -= 1
        >>> all(
        ...   HR(n).rank
        ...   in admissible_lorentz_ranks(HR(n).dim)
        ...   for n in range(n_max+1)
        ... )
        True

    In fact, we know the formula for at least one such simulacrum::

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
    # cached cones (the ``simulacra`` method uses it implicitly).
    n_max = 0
    mcd = max_cone_dim()
    while HR(n_max).dim <= mcd:
        n_max += 1
    n_max -= 1

    n_min = 3
    if n_max <= n_min:
        # No cached cones?
        return True

    return all( HR(n).simulacra() for n in range(n_min, n_max+1) )



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

    From the proof of Proposition 5, we know that for ``n >= 4``,
    ``HC(n)`` has symmetric simulacra with only Lorentz
    factors. Moreover when ``n < 3``, ``HC(n)`` _is_ a Lorentz
    cone. In either case, ``HC(n)`` should share its signature with a
    sum of Lorentz cones. We begin by computing the largest ``n`` for
    which we have the corresponding sum-of-Lorentz-cone data cached
    (it's easy in this case)::

        >>> from math import floor, sqrt
        >>> from sql import admissible_lorentz_ranks, max_lorentz_rank_dim
        >>> n_max = -1
        >>> mlrd = max_lorentz_rank_dim()
        >>> if mlrd >= 0:
        ...     n_max = floor(sqrt(mlrd))
        >>> all(
        ...   HC(n).rank
        ...   in admissible_lorentz_ranks(HC(n).dim)
        ...   for n in range(n_max+1)
        ...   if n != 3
        ... )
        True

    We know the simulacra for ``n >= 4`` explicitly; they are given in
    the proof of the proposition::

        >>> all(
        ...   HC(n).signature() == K.signature()
        ...   for n in range(4,100)
        ...   if (K2 := RN(n**2 - 5*n + 6))
        ...   and (K := DirectSum(2*[L(n+1)] + [K2, L(4)] + (n-4)*[L(3)]))
        ... )
        True

    An example using the :meth:`SymmetricCone.simulacra` method for
    ``n = 4``::

        >>> from sql import max_cone_dim
        >>> K = DirectSum([L(5), L(5), L(4), RN(2)])
        >>> mcd = max_cone_dim()
        >>> skip = ( HC(4).dim > mcd )
        >>> skip or K in HC(4).simulacra()
        True

    """
    from math import floor, sqrt
    from sql import max_cone_dim

    # Figure out how big "n" can be if we want to use the database of
    # cached cones (the ``simulacra`` method uses it implicitly).
    n_max = 0
    mcd = max_cone_dim()
    if mcd >= 0:
        n_max = floor(sqrt(mcd))

    n_min = 4
    if n_max <= n_min:
        # No cached cones?
        return True

    return ( not HC(3).simulacra()
             and
             all(HC(n).simulacra() for n in range(n_min, n_max+1)) )



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

    From the proof of Proposition 6, we know that for ``n >= 3``,
    ``HH(n)`` has symmetric simulacra with only Lorentz
    factors. Moreover when ``n < 3``, ``HH(n)`` _is_ a Lorentz
    cone. In either case, ``HH(n)`` should share its signature with a
    sum of Lorentz cones. We begin by computing the largest ``n`` for
    which we have the corresponding sum-of-Lorentz-cone data cached::

        >>> from sql import admissible_lorentz_ranks, max_lorentz_rank_dim
        >>> n_max = -1
        >>> mlrd = max_lorentz_rank_dim()
        >>> while HH(n_max).dim <= mlrd:
        ...     n_max += 1
        >>> n_max -= 1
        >>>
        >>> all(
        ...   HH(n).rank
        ...   in admissible_lorentz_ranks(HH(n).dim)
        ...   for n in range(n_max+1)
        ... )
        True

    We know the formula for simulacra explicitly; they are given in
    the proof of the proposition::

        >>> all(
        ...   HH(n).signature() == K.signature()
        ...   for n in range(3,100)
        ...   if (K3 := RN(2*n**2 - 3*n - 2))
        ...   and (K := DirectSum([L(2*n+2), K3]))
        ... )
        True

    Further checks of the low-dimensional formulas using the
    :meth:`SymmetricCone.simulacra` method::

        >>> from sql import max_cone_dim
        >>> mcd = max_cone_dim()

        >>> K = DirectSum([L(8), RN(7)])
        >>> skip = ( HH(3).dim > mcd )
        >>> skip or K in HH(3).simulacra()
        True

        >>> K = DirectSum([L(10), RN(18)])
        >>> skip = ( HH(4).dim > mcd )
        >>> skip or K in HH(4).simulacra()
        True

        >>> K = DirectSum([L(12), RN(33)])
        >>> skip = ( HH(5).dim > mcd )
        >>> skip or K in HH(5).simulacra()
        True

    """
    from sql import max_cone_dim

    # Figure out how big "n" can be if we want to use the database of
    # cached cones (the ``simulacra`` method uses it implicitly).
    n_max = 0
    mcd = max_cone_dim()
    while HH(n_max).dim <= mcd:
        n_max += 1
    n_max -= 1

    n_min = 3
    if n_max <= n_min:
        # No cached cones?
        return True

    return all( HH(n).simulacra() for n in range(n_min, n_max+1) )


def test_proposition7() -> bool:
    r"""
    Test the statement of Proposition 7.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The result should hold::

        >>> test_proposition7()
        True

    From the proof of Proposition 7, we know that ``HO(3)`` has a
    symmetric simulacrum with only Lorentz factors. Moreover when ``n
    < 3``, ``HO(n)`` _is_ a Lorentz cone. In either case, ``HO(n)``
    should share its signature with a sum of Lorentz cones. We begin
    by computing the largest ``n`` for which we have the corresponding
    sum-of-Lorentz-cone data cached::

        >>> from sql import admissible_lorentz_ranks, max_lorentz_rank_dim
        >>> n_max = -1
        >>> mlrd = max_lorentz_rank_dim()
        >>> while n_max <= 3 and HO(n_max).dim <= mlrd:
        ...     n_max += 1
        >>> n_max -= 1
        >>>
        >>> all(
        ...   HO(n).rank
        ...   in admissible_lorentz_ranks(HO(n).dim)
        ...   for n in range(n_max+1)
        ... )
        True

    We know the formula for one simulacrum explicitly; it is given
    in the proof of the proposition::

        >>> K = DirectSum([L(11),L(5),L(3),RN(8)])
        >>> K.signature() == HO(3).signature()
        True

    Repeat with the cached simulacra data::

        >>> from sql import have_cone_dim
        >>> (not have_cone_dim(HO(3).dim)) or K in HO(3).simulacra()
        True

    """
    from sql import max_cone_dim

    if HO(3).dim > max_cone_dim():
        # no data
        return True

    return not (not HO(3).simulacra())


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
    return all ( not (not K.simulacra())
                 or K == HC(3)
                 or isinstance(K,L)
                 for K in Ks
                 if K.dim <= max_cone_dim() )


def test_lemma2() -> bool:
    r"""
    Test the statement of Lemma 2.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The result should hold::

        >>> test_lemma2()
        True

    """
    from sql import max_cone_dim
    Ks = ( random_cone() for _ in range(100) )

    return all(
      K.simulacra()
      or
      all(not K_i.simulacra() for K_i in K.factors())
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
      not (not K.simulacra())
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
    r |= any( partition_rank(p) == K.rank
              for p in partitions(K.dim, include_two=False) )

    dim_rest = K.dim - 9
    if dim_rest >= 0:
        # Try with an HC(3) factor
        r |= any( partition_rank(p) == (K.rank - 17)
                  for p in partitions(dim_rest, include_two=False) )
    return r


def test_theorem3_example() -> bool:
    r"""
    Test the example given directly after Theorem 3.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The relationships in the example should hold::

        >>> test_theorem3_example()
        True

    """
    result = True

    # The two simulacra of HC(3) mentioned in the example
    K1 = DirectSum([L(11),L(3)] + [L(5),RN(8)])
    K2 = DirectSum([L(11),L(3)] + [L(4)] + 3*[L(3)])

    # Signature comparison, suffices because they're obviously
    # not isomorphic
    result &= ( K1.signature() == K2.signature() )
    result &= ( K1.signature() == HO(3).signature() )

    # And repeat using cached simulacra if possible
    from sql import have_cone_dim
    if have_cone_dim(HO(3).dim):
        result &= K1 in K2.simulacra()
        result &= K2 in K1.simulacra()

        result &= K1 in HO(3).simulacra()
        result &= HO(3) in K1.simulacra()

        result &= K2 in HO(3).simulacra()
        result &= HO(3) in K2.simulacra()

    return result


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

    Derive the result by minimizing ``K.rank`` (for a fixed ``K.dim``)
    in :func:`lowerbound2`::

        >>> from sympy import symbols
        >>> n,d = symbols("n,d", integer=True, positive=True)
        >>> fix_floor(n >= 2 + L(d+1).rank - L(d).rank).expand()
        n >= d + 2

    """
    from random import randint
    K_n_pairs = ( (random_cone(),randint(0,1000))
                  for _ in range(100) )
    return all(
      (n >= K.dim + 2) or any([n < lowerbound1(K),
                               n < lowerbound2(K),
                               n < lowerbound3(K)])
      for (K,n) in K_n_pairs
    )


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
        Ln_plus_K = DirectSum([K, L(n)])

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



def test_lemma5() -> bool:
    r"""
    Test Lemma 5.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The implication in the result should hold::

        >>> test_lemma5()
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
        >>> from partitions import f
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

    Ks = ( random_cone() for _ in range(10) )

    result = True
    for K in Ks:
        # Partitioning 60 can be done in a few seconds, but e.g. 80
        # may crash the machine. Keep in mind that we're partitioning
        # n + dim(K), and that this is multiplied (at worst) by the
        # length of Ks!
        min_n = max(lowerbound1(K), lowerbound2(K), lowerbound3(K))
        max_n = 60 - K.dim

        if min_n > max_n:
            # This is a pretty tight window. For example we know that
            # min_n = 31 for K == HC(3), but then K.dim == 9 already.
            continue
        n = randint(min_n, max_n)

        J = DirectSum([K,L(n)])
        for p in partitions(n + K.dim, n-1, include_two=False):
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
    irreducible factors with sums of Lorentz cones. The replacements
    should agree in dimension, and have a dominant Lyapunov rank. The
    first example we give is for the 3x3 complex PSD cone::

        >>> L(9).dim == HC(3).dim
        True
        >>> L(9).rank > HC(3).rank
        True

    For all other non-Lorentz factors, we verify explicitly that any
    of their all-Lorentz simulacra will have no factor of dimension
    ``n`` or greater::

        >>> from sql import max_cone_dim
        >>> K = random_cone()
        >>> min_n = max(lowerbound1(K), lowerbound2(K), lowerbound3(K))
        >>> max_n = max_cone_dim() - K.dim
        >>> all(
        ...   (not isinstance(f,L)) or f.dim < n
        ...   for n in range(min_n, max_n)
        ...   for J in DirectSum([K,L(n)]).simulacra()
        ...   for I in J.factors()
        ...   if (not isinstance(I,L)) and I.dim >= n
        ...   for I_prime in I.simulacra()
        ...   for f in I_prime.factors()
        ... )
        True

    As part of the argument, we claim that all real, complex, and
    quaternion PSD factors ``I`` satisfy the bound ``12*I.dim >=
    5*I.rank``. We confirm this with a random sample of irreducible
    cones, and set ``n`` large enough to ensure that the factors we
    generate are not isomorphic to Lorentz cones::

        >>> all(
        ...   12*K.dim >= 5*K.rank
        ...   for _ in range(100)
        ...   if (K := random_irreducible_cone(3,20))
        ...   and not isinstance(K, (L,HO))
        ... )
        True

    ...Although, this is easier to just check::

        >>> from sympy import symbols
        >>> n = symbols("n", integer=True, positive=True)
        >>> K = SymmetricCone(0)
        >>> # HR(n)
        >>> K.dim = n**2
        >>> K.dim = (n**2 + n)/2
        >>> K.rank = n**2
        >>> 12*K.dim - 5*K.rank
        n**2 + 6*n
        >>> # HC(n)
        >>> K.dim = n**2
        >>> K.rank = 2*n**2 - 1
        >>> 12*K.dim - 5*K.rank
        2*n**2 + 5
        >>> # HH(n)
        >>> K.dim = 2*n**2 - n
        >>> K.rank = 4*n**2
        >>> 12*K.dim - 5*K.rank
        4*n**2 - 12*n

    We also claim that ``5*L(n).rank >= 12*(2*n - 1)`` as part of the
    argument.  This holds under our assumption that ``n >= 10``, and
    no smaller bound will work::

        >>> all(
        ...   5*L(n).rank >= 12*(2*n - 1)
        ...   for n in range(10, 100)
        ... )
        True
        >>> 5*L(9).rank >= 12*(2*9 - 1)
        False

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

    Test the claim that if we drop the ``n >= 10`` condition on
    Theorem 4, then the only new counterexample we get is ``HR(3) ~
    L(4) + L(2)`` at ``n == 4``. Thus the real condition is ``n !=
    4``, for which we get no counterexamples::

        >>> from sql import all_cones_of_dim, max_cone_dim
        >>> actual = []
        >>> # skip check if not enough data; the largest dimension
        >>> # we need is 14 when n=9 and d=5.
        >>> skip = ( max_cone_dim() < 14 )
        >>> if not skip:
        ...     for n in range(10):
        ...         for d in range(max_dimK(n)+1):
        ...             for K in all_cones_of_dim(d):
        ...                 if n < lowerbound1(K) or n < lowerbound2(K):
        ...                     continue
        ...                 C = DirectSum([K, L(n)])
        ...                 for s in C.simulacra():
        ...                     if L(n) not in s.factors():
        ...                         actual.append( (n, (C,s)) )
        >>>
        >>> expected = [ (4, (DirectSum([L(2),L(4)]), HR(3))) ]
        >>> skip or (actual == expected)
        True

    For a given ``K.dim``, :func:`lowerbound1` is minimized by the
    nonnegative orthant (with value ``n >= 2``), and
    :func:`lowerbound2` is minimized by the Lorentz cone with value
    ``K.dim + 2`` (proof: we are either maximizing beta(K) or
    minimizing the Lyapunov rank in fixed dimensions). The second
    bound is therefore increasing with ``K.dim`` and obviously
    dominates the first. As a result, we can use :func:`lowerbound2`
    to determine the first potentially valid ``n`` corresponding to
    any ``K.dim``. Moreover we can use this to compute the first
    ``K.dim`` such that ``n + K.dim`` will exceed the largest
    dimension we have cached: basically we just add ``K.dim`` to the
    bound and set it greater than the largest dimension we have
    cached, i.e. we solve ``K.dim + 2 + K.dim > max_cone_dim()``::

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

        >>> from itertools import chain
        >>> from sympy import symbols
        >>> n,d = symbols("n,d", integer=True, positive=True)
        >>> f = lambda x: (x**2 - x + 2)/2   # no integer division... sympy
        >>> def rank_too_small(X):
        ...     g = f(n) + d - X.rank - f(n+d-X.dim)
        ...     return all( g.subs({n:i,d:j}) > 0
        ...                 for i in chain(range(0,4), range(5,10))
        ...                 for j in range(max_dimK(i)+1)
        ...                 if i+j-X.dim >= 0 )

    The argument we use to rule out ``HR(3)``, ``HR(4)``, and
    ``HC(3)`` factors::

        >>> X = HR(3)
        >>> rank_too_small(X)
        True

        >>> X = HR(4)
        >>> rank_too_small(X)
        True

        >>> X = HC(3)
        >>> rank_too_small(X)
        True

    One of the last statements in the proof is that the conclusion is
    easy to verify for the "new" cases because there simply aren't any
    new simulacra (so we don't even have to worry about whether or not
    "J" has the stated form). We're checking that a list of simulacra
    is empty, so there is no need to check the dimension of our cone
    against :func:`sql.max_cone_dim`; if we exceed it, we'll get back
    empty lists of simulacra anyway::

        >>> from sql import all_cones_of_dim
        >>> all( not DirectSum([K,L(n)]).simulacra()
        ...      for n in chain(range(0,4), range(5,10))
        ...      for d in range(max_dimK(n)+1)
        ...      for K in all_cones_of_dim(d)
        ...      if  n >= max(lowerbound1(K),lowerbound2(K)) )
        True

    Finally, we check the proof using the low-tech method that we have
    described: partitions. There are many matching signatures, but
    they're all from isomorphic cones once you consider that ``L(2) ==
    RN(2)``::

        >>> from partitions import partitions, f
        >>>
        >>> # We'll collect the matching signatures in a list
        >>> matches = []
        >>>
        >>> # reimplement the lower bounds in terms of partitions
        >>> lb1 = lambda p: 2 + partition_rank(p) - sum(p)
        >>> lb2 = lambda p: 2 + f(1+sum(p)) - partition_rank(p)
        >>>
        >>> for n in chain(range(0,4), range(5,10)):
        ...     # start at d=1 to avoid getting [0,n] ~ [n]
        ...     for d in range(1, max_dimK(n)+1):
        ...         for K in partitions(d, include_two=False):
        ...             if n < lb1(K): continue
        ...             if n < lb2(K): continue
        ...             Ln_K_rank = f(n) + partition_rank(K)
        ...             for J in partitions(d+n, include_two=False):
        ...                 if (partition_rank(J) == Ln_K_rank):
        ...                         # n is guaranteed to be larger than d,
        ...                         # so we know it goes at the end. Special
        ...                         # case to avoid inserting zero at the
        ...                         # beginning of a partition.
        ...                         if K == [0]:
        ...                             Ln_K = [n]
        ...                         else:
        ...                             Ln_K = K + [n]
        ...                         if not J == Ln_K:
        ...                             matches.append( (Ln_K, J) )
        >>> matches
        []

    Check the relationships between the lower bounds on ``n``. Each
    can be violated while the others are satisfied, and in each case,
    Theorem 5 fails. This is discussed subsequent to Theorem 5 in the
    paper::

        >>> from random import randint
        >>> m = randint(5,20)
        >>> K = L(m)
        >>> lowerbound2(K) < lowerbound1(K)
        True
        >>> n = lowerbound1(K) - 1
        >>> (n >= lowerbound1(K), n >= lowerbound2(K), n != 4)
        (False, True, True)
        >>> ( DirectSum([L(m), L(n)]).signature()
        ...   ==
        ...   DirectSum([RN(m-1), L(n+1)]).signature() )
        True

        >>> K = RN(2)
        >>> n = 4
        >>> lowerbound1(K) < lowerbound2(K)
        True
        >>> (n >= lowerbound1(K), n >= lowerbound2(K), n != 4)
        (True, True, False)
        >>> DirectSum([K,L(n)]).signature() == HR(3).signature()
        True

        >>> K = HC(3)
        >>> lowerbound1(K) < lowerbound2(K)
        True
        >>> n = 30
        >>> (n >= lowerbound1(K), n >= lowerbound2(K), n != 4)
        (True, False, True)
        >>> ( DirectSum([K,L(n)]).signature()
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
    # will appear to have no simulacra, and (for the sake of the
    # theorem) that's fine. But we do have to stop _somewhere_, so we
    # might as well stop here?
    max_d = (max_cone_dim() - 2) // 2

    for d in range(1, max_d+1):
        min_n = max(lowerbound1(RN(d)), lowerbound2(L(d)))
        max_n = max_cone_dim() - d
        for n in range(min_n, max_n+1):
            if n == 4: continue
            for K in all_cones_of_dim(d):
                if n < lowerbound1(K): continue
                if n < lowerbound2(K): continue
                lhs = DirectSum([K,L(n)])
                for J in lhs.simulacra():
                    fs = list(J.factors())
                    # remove() raises an error if L(n) isn't a factor,
                    # so this guarantees that L(n) is one.
                    fs.remove(L(n))
                    J_prime = DirectSum(fs)
                    result &= J_prime in K.simulacra()

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
        >>> n_max = 0
        >>> n_max = max_cone_dim() - 9
        >>> all( not DirectSum([HC(3),L(n)]).simulacra()
        ...      for n in range(n_min, n_max+1) )
        True

    If you think hard about it, or consult an earlier version of the
    paper, you will conclude that only sums of Lorentz cones need to
    be checked for simulacra of ``HC(3) + L(n)``. Here we repeat the
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
            min_n = max(lowerbound1(K),lowerbound2(K))
            if min_n > max_n:
                # The dimension of L(n)+K is guaranteed to exceed the
                # dimensions we have cached, so skip this cone.
                continue
            if min_n == 4 and max_n == 4:
                # We're gonne get stuck generating random 4s (because
                # 4 is not a valid value for n) forever
                continue
            if not K.simulacra():
                n = 4
                while n == 4:
                    n = randint(min_n, max_n)
                # the corollary says that L(n)+K should have no
                # simulacra
                result &= not DirectSum([K,L(n)]).simulacra()

    return result


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

    # the list of "n" where we expect to find simulacra
    expected_n = [2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 18, 21, 22, 30]

    # We're going to be calling simulacra() on L(n)+HC(3), so make
    # sure we don't make "n" so large that we exceed what is cached.
    n_max = max_cone_dim() - 9

    return all(
      (not DirectSum([HC(3),L(n)]).simulacra())
      or
      (n in expected_n)
      for n in range(n_max+1)
    )


def test_proposition9() -> bool:
    r"""
    Test Proposition 9.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The implication in the result should hold::

        >>> test_proposition9()
        True

    Find the pairs of ``m`` and ``n`` that aren't a corollary of
    Theorem 5. The cases where ``m == 0`` are trivial::

        >>> [
        ...   (m,n)
        ...   for m in range(1,100)
        ...   for n in range(100)
        ...   if (
        ...     n < lowerbound1(L(m))
        ...     or n < lowerbound2(L(m))
        ...     or n == 4
        ...   )
        ...   and n >= (m**2 - 3*m + 6)/2
        ... ]
        [(1, 2), (1, 4), (2, 2), (2, 3), (2, 4), (3, 3), (3, 4), (4, 5)]

    Verify the cases mentioned explicitly in the proof. First, the
    ``m != 2`` cases where there are no simulacra::

        >>> from sql import max_cone_dim
        >>> mcd = max_cone_dim()

        >>> m = 4
        >>> K = L(m)
        >>> n = 5
        >>> (n >= lowerbound1(K), n >= lowerbound2(K), n != 4)
        (True, False, True)
        >>> if K.dim + n <= mcd:
        ...     DirectSum([K,L(n)]).simulacra()
        ... else:
        ...     ()
        ()

        >>> m = 3
        >>> K = L(m)
        >>> n = 4
        >>> (n >= lowerbound1(K), n >= lowerbound2(K), n != 4)
        (True, False, False)
        >>> if K.dim + n <= mcd:
        ...     DirectSum([K,L(n)]).simulacra()
        ... else:
        ...     ()
        ()
        >>> n = 3
        >>> (n >= lowerbound1(K), n >= lowerbound2(K), n != 4)
        (True, False, True)
        >>> if K.dim + n <= mcd:
        ...     DirectSum([K,L(n)]).simulacra()
        ... else:
        ...     ()
        ()

        >>> m = 1
        >>> K = L(m)
        >>> n = 2
        >>> (n >= lowerbound1(K), n >= lowerbound2(K), n != 4)
        (True, False, True)
        >>> if K.dim + n <= mcd:
        ...     DirectSum([K,L(n)]).simulacra()
        ... else:
        ...     ()
        ()
        >>> n = 4
        >>> (n >= lowerbound1(K), n >= lowerbound2(K), n != 4)
        (True, True, False)
        >>> if K.dim + n <= mcd:
        ...     DirectSum([K,L(n)]).simulacra()
        ... else:
        ...     ()
        ()

    And now the ``m == 2`` case where there is exactly one
    counterexample::

        >>> from sql import max_cone_dim
        >>> mcd = max_cone_dim()
        >>> m = 2
        >>> K = L(m)
        >>> ns = [2, 3]
        >>> any( n < lowerbound1(K) for n in ns )  # not violated
        False
        >>> all( n < lowerbound2(K) for n in ns )  # violated
        True
        >>> sims = []
        >>> for n in ns:
        ...     if K.dim + n <= mcd:
        ...         sims.append(DirectSum([K,L(n)]).simulacra())
        ...     else:
        ...         sims.append(())
        >>> sims
        [(), ()]

        >>> from sql import max_cone_dim
        >>> expected = (HR(3),)
        >>> J = DirectSum([K,L(4)])
        >>> skip = ( J.dim > max_cone_dim() )
        >>> skip or (J.simulacra() == expected)
        True

   Confirm that the stated bound actually comes from
   :func:`lowerbound1`::

        >>> from sympy import symbols
        >>> m = symbols("m", integer=True, positive=True)
        >>> fix_floor(lowerbound1(L(m))) == (m**2 - 3*m + 6)/2
        True

    The ``L(m) + L(n)`` cases with ``m,n <= 2`` are singled out in
    a bullet point::

        >>> all(
        ...   not list(partition_simulacra(p))
        ...   for k in range(5)
        ...   for p in partitions(k)
        ... )
        True

    The ``L(3) + L(2)`` case is also singled out in a bullet point::

        >>> list(partition_simulacra([3,2]))
        []

    """
    from sql import max_cone_dim
    mcd = max_cone_dim()

    return all(
      not DirectSum([L(m),L(n)]).simulacra()
      for m in range(1, mcd + 1)
      for n in range(lowerbound1(L(m)), mcd - m)
      if m != 2
    )


def test_lemma6() -> bool:
    r"""
    Test Lemma 6.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The implication in the result should hold::

        >>> test_lemma6()
        True

    For ``n <= 2``, there shouldn't be any simulacra in the first
    place::

        >>> from sql import max_cone_dim
        >>> mcd = max_cone_dim()
        >>> [ K.simulacra() if 2*n <= mcd else ()
        ...   for n in range(3)
        ...   if (K := DirectSum(2*[L(n)])) ]
        [(), (), ()]

    """
    from sql import all_cones_of_dim, max_cone_dim

    # For ALL cones with simulacra, there EXISTS a simulacrum, such
    # that ALL of its factors are HC(3) or Lorentz cones.
    n_max = max_cone_dim() // 2
    return all(
      not K.simulacra()
      or
      any(
        all( f == HC(3) or isinstance(f,L) for f in J.factors() )
        for J in K.simulacra()
      )
      for n in range(n_max+1)
      if (K := DirectSum(2*[L(n)]))
    )


def test_proposition10() -> bool:
    r"""
    Test Proposition 10.

    Returns
    -------

    ``True`` if the test passed, and ``False`` otherwise.

    Examples
    --------

    The characterization in the result should hold::

        >>> test_proposition10()
        True

    In addition to the ``simulacra`` check, we can also use our cached
    partitions thanks to Lemma 6. This allows us to test all the way
    up to ``n == 100``.....


    Define the "constants" we use in the proof::

        >>> m = lambda n: (n // 5)
        >>> k = lambda n: n - 5*m(n)
        >>> r = lambda n: m(n) - k(n)**2 + 1 - 3*((m(n) - k(n)**2 + 1) // 3)
        >>> alpha = lambda n: (m(n) - 4*k(n)**2 + 15*k(n) - 14 - r(n)) / 3
        >>> gamma = lambda n: 2*m(n) - 22*k(n) - (4*m(n) - 16*k(n)**2 - 68 + 5*r(n))/3

    And check some of their properties::

        >>> n_max = 200
        >>> all( n == 5*m(n) + k(n) for n in range(n_max+1) )
        True
        >>> all( alpha(n).is_integer() for n in range(n_max+1) )
        True
        >>> all( gamma(n).is_integer() for n in range(n_max+1) )
        True
        >>> all( alpha(n) >= 0 for n in range(100, n_max+1) )
        True
        >>> all( gamma(n) >= 0 for n in range(100, n_max+1) )
        True

    Finally, we check the dimension and Lyapunov rank of our
    simulacrum symbolically::

        >>> from sympy import expand, symbols
        >>> from partitions import f
        >>> n,m = symbols("n,m", integer=True, positive=True)
        >>> k,r = symbols("k,r", integer=True, nonnegative=True)
        >>> alpha = (m - 4*k**2 + 15*k - 14 - r) / 3
        >>> gamma = 2*m - 22*k - (4*m - 16*k**2 - 68 + 5*r)/3
        >>> J_dim = (7*m + k) + (m + 3*k - 4) + 4*alpha + 3*r + gamma
        >>> J_rank = fix_floor(f(7*m + k)) + fix_floor(f(m + 3*k - 4)) + alpha*f(4) + r*f(3) + gamma
        >>> K_dim = 2*n
        >>> K_rank = 2*fix_floor(f(n))
        >>> (K_dim - J_dim).subs({n: 5*m + k})
        0
        >>> expand((K_rank - J_rank).subs({n: 5*m + k}))
        0

    Since the target cone has exactly two factors, it suffices (per
    the proof) to check partitions with/without an ``HC(3)`` offset.
    We do this only up to ``n == 18`` because all greater ``n`` lead
    to simulacra, and the existence of a simulacra is much easier to
    verify by just writing down its factors::

        >>> from partitions import partition_rank, partition_simulacra
        >>> n_max = 18
        >>> n_without_simulacra = []
        >>> for n in range(n_max+1):
        ...     p = [n,n]
        ...     simcount = len(list(partition_simulacra(p)))
        ...     f = lambda q: partition_rank(q) == (partition_rank(p) - 17)
        ...     if n >= 5:
        ...         simcount += len(list(
        ...           filter(f, partitions(2*n - 9, include_two=False))
        ...         ))
        ...     if simcount == 0:
        ...         n_without_simulacra.append(n)
        >>> n_without_simulacra
        [0, 1, 2, 3, 5, 6, 7, 11, 12, 13, 18]

    Finally, we demonstrate simulacra for all ``n`` between ``0`` and
    ``100`` that are not in the no-simulacra list::

        >>> d = {
        ...   4: [1, 2, 5],
        ...   8: [1, 5, 10],
        ...   9: [1, 2, 3, 12],
        ...   10: [2, 5, 13],
        ...   14: [1, 10, 17],
        ...   15: [1, 3, 6, 20],
        ...   16: [2, 2, 2, 2, 2, 22],
        ...   17: [2, 4, 5, 23],
        ...   19: [1, 2, 2, 2, 5, 26],
        ...   20: [2, 13, 25],
        ...   21: [1, 2, 12, 27],
        ...   22: [1, 17, 26],
        ...   23: [1, 3, 12, 30],
        ...   24: [1, 2, 2, 2, 2, 6, 33],
        ...   25: [2, 3, 12, 33],
        ...   26: [5, 13, 34],
        ...   27: [2, 2, 2, 12, 36],
        ...   28: [1, 5, 13, 37],
        ...   29: [1, 2, 2, 2, 2, 2, 7, 40],
        ...   30: [1, 2, 3, 4, 9, 41],
        ...   31: [3, 8, 9, 42],
        ...   32: [1, 26, 37],
        ...   33: [1, 3, 20, 42],
        ...   34: [2, 25, 41],
        ...   35: [2, 8, 13, 47],
        ...   36: [3, 3, 19, 47],
        ...   37: [1, 2, 2, 5, 14, 50],
        ...   38: [1, 2, 4, 19, 50],
        ...   39: [1, 2, 27, 48],
        ...   40: [2, 2, 4, 19, 53],
        ...   41: [2, 4, 23, 53],
        ...   42: [1, 2, 3, 3, 19, 56],
        ...   43: [2, 5, 23, 56],
        ...   44: [1, 37, 50],
        ...   45: [1, 3, 30, 56],
        ...   46: [5, 29, 58],
        ...   47: [1, 2, 2, 2, 26, 61],
        ...   48: [1, 2, 2, 2, 6, 18, 65],
        ...   49: [2, 2, 4, 26, 64],
        ...   50: [1, 2, 10, 20, 67],
        ...   51: [2, 3, 33, 64],
        ...   52: [2, 41, 61],
        ...   53: [2, 5, 31, 68],
        ...   54: [1, 3, 4, 30, 70],
        ...   55: [2, 2, 2, 5, 26, 73],
        ...   56: [10, 29, 73],
        ...   57: [2, 2, 2, 36, 72],
        ...   58: [1, 50, 65],
        ...   59: [1, 3, 42, 72],
        ...   60: [1, 2, 2, 2, 2, 33, 78],
        ...   61: [4, 5, 34, 79],
        ...   62: [1, 2, 3, 4, 33, 81],
        ...   63: [1, 2, 48, 75],
        ...   64: [1, 2, 2, 44, 79],
        ...   65: [4, 7, 34, 85],
        ...   66: [1, 2, 4, 5, 33, 87],
        ...   67: [3, 3, 4, 37, 87],
        ...   68: [1, 2, 7, 38, 88],
        ...   69: [1, 2, 2, 6, 37, 90],
        ...   70: [3, 3, 47, 87],
        ...   71: [2, 4, 8, 34, 94],
        ...   72: [1, 2, 3, 4, 41, 93],
        ...   73: [1, 2, 2, 2, 2, 2, 40, 95],
        ...   74: [1, 65, 82],
        ...   75: [1, 3, 56, 90],
        ...   76: [1, 2, 4, 50, 95],
        ...   77: [2, 4, 53, 95],
        ...   78: [1, 3, 6, 46, 100],
        ...   79: [2, 3, 57, 96],
        ...   80: [5, 58, 97],
        ...   81: [1, 4, 7, 45, 105],
        ...   82: [2, 2, 4, 53, 103],
        ...   83: [2, 5, 56, 103],
        ...   84: [3, 3, 59, 103],
        ...   85: [4, 7, 50, 109],
        ...   86: [10, 53, 109],
        ...   87: [2, 3, 64, 105],
        ...   88: [5, 65, 106],
        ...   89: [1, 2, 2, 2, 61, 110],
        ...   90: [1, 4, 6, 54, 115],
        ...   91: [1, 7, 9, 45, 120],
        ...   92: [1, 82, 101],
        ...   93: [1, 2, 75, 108],
        ...   94: [1, 2, 2, 2, 2, 61, 118],
        ...   95: [2, 2, 4, 64, 118],
        ...   96: [1, 2, 4, 5, 57, 123],
        ...   97: [2, 5, 68, 119],
        ...   98: [2, 3, 8, 57, 126],
        ...   99: [2, 2, 2, 72, 120],
        ...   100: [2, 85, 113]
        ... }
        >>> all(
        ...   K.rank == partition_rank(d[n])
        ...   for n in range(101)
        ...   if not n in n_without_simulacra
        ...   if (K := DirectSum(2*[L(n)]))
        ... )
        True

    """
    from sql import max_cone_dim
    n_max = max_cone_dim() // 2
    n_without_simulacra = [
      n for n in range(n_max+1)
      if (K := DirectSum([L(n)]*2))
      and not K.simulacra()
    ]

    # We have to filter the expected result based on the
    # max available cone dimension, too.
    expected = [ e for e in [0, 1, 2, 3, 5, 6, 7, 11, 12, 13, 18]
                 if e <= n_max ]

    return (n_without_simulacra == expected)
