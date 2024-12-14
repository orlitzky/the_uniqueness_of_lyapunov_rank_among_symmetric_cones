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

Check the symbolic identity in Proposition 8 for the derivative of the
g-delta function::

    >>> from sympy import diff, symbols
    >>> x,d = symbols("x,d", integer=True, positive=True)
    >>> g = fix_floor(L(x).rank - L(x-d).rank)
    >>> diff(g, x)
    d

Check the symbolic identity in Proposition 8 for the Lyapunov rank of
``L(n-1)`` in terms of that of ``L(n)``::

    >>> from sympy import simplify, symbols
    >>> n = symbols("n", integer=True, positive=True)
    >>> lhs = L(n-1).rank
    >>> rhs = L(n).rank - (n-1)
    >>> simplify(fix_floor(lhs - rhs))
    0

Check implication (3) in Proposition 8::

    >>> f = lambda x: L(x).rank
    >>> results = []
    >>> for n in range(100):
    ...     for d in range(100):
    ...         r = (n-1)<(10+d) or f(n-1)+f(10) >= f(n-1-d)+f(10+d)
    ...         results.append(r)
    >>> all(results)
    True

Check the relationships between the lower bounds on ``n``::

    >>> K = L(3)
    >>> lowerbound1(K) < lowerbound2(K) < lowerbound3(K)
    True

    >>> K = RN(5)
    >>> lowerbound1(K) < lowerbound3(K) < lowerbound2(K)
    True

    >>> K = L(6)
    >>> lowerbound2(K) < lowerbound3(K) < lowerbound1(K)
    True


Test the ``n < 9`` case of Giovanni's Theorem 3::

    >>> def check(m,n):
    ...     d = m+n
    ...     if d < 11:
    ...         # no similacra, RHS is too big, and we need at least
    ...         # two lorentz factors to avoid the previous result
    ...         return True
    ...     ranks = ( 17 + r for r in admissible_lorentz_ranks(d - 9) )
    ...     return (f(m)+f(n)) not in ranks
    ...
    >>> all( check(m,n) for m in range(1,9) for n in range(max(2,m),9) )
    True

Check Giovanni's Theorem 3 directly::

    >>> def check(m,n):
    ...     K = DirectSum([L(m),L(n)])
    ...     s = K.similacra()
    ...     if not s:
    ...         return True
    ...     if len(s) == 1 and s[0] == HR(3):
    ...         return True
    ...     ub = (m**2 - 3*m + 4) // 2
    ...     return ub >= n >= m >= 4
    >>> from sql import max_cone_dim
    >>> # need m+n <= max_cone_dim()
    >>> m_max = max_cone_dim() // 2
    >>> all( check(m,n)
    ...      for m in range(4, m_max)
    ...      for n in range(m, max_cone_dim() - m) )
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

Check Theorem 4 using our precomputed dictionary of cones. We start
with a cone ``K``, and then add ``L(n)`` factors to it. If the
resulting sum has similacra, then each similacrum should have an
``L(n)`` factor. We do this for as many ``n`` as we can, constrained
by the fact that ``K + L(n)`` needs to have a dimension that we have
cached::

    >>> from sql import max_cone_dim
    >>> def check(K):
    ...     n_min = max(lowerbound1(K),lowerbound2(K),lowerbound3(K))
    ...     n_max = max_cone_dim() - K.dim
    ...     result = True
    ...     for n in range(n_min, n_max+1):
    ...         lhs = DirectSum([L(n),K])
    ...         result &= all( L(n) in f.factors()
    ...                        for f in lhs.similacra() )
    ...     return result

Some small, hand-crafted examples::

   >>> K = DirectSum([HC(3), L(1)])
   >>> check(K)
   True

   >>> K = DirectSum([HC(3), L(2)])
   >>> check(K)
   True

   >>> K = DirectSum([HC(3), L(3)])
   >>> check(K)
   True

   >>> K = DirectSum([HC(3), L(4)])
   >>> check(K)
   True

   >>> K = DirectSum([HC(3), L(3), L(1)])
   >>> check(K)
   True

   >>> K = DirectSum([HC(3), RN(4)])
   >>> check(K)
   True

   >>> K = DirectSum([HR(3), RN(4)])
   >>> check(K)
   True

   >>> K = DirectSum([HR(3), L(1)])
   >>> check(K)
   True

   >>> K = DirectSum([HR(3), L(2)])
   >>> check(K)
   True

   >>> K = DirectSum([HR(3), L(3)])
   >>> check(K)
   True

   >>> K = DirectSum([HR(3), HR(3)])
   >>> check(K)
   True

   >>> K = DirectSum([L(3)]*2)
   >>> check(K)
   True

   >>> K = DirectSum([L(3)]*3)
   >>> check(K)
   True

   >>> K = RN(5)
   >>> check(K)
   True

A random example. Unfortunately the lower bound on ``n`` is almost
always going to be higher than the largest dimension we have cached,
and this will be a no-op in that case::

    >>> from cones import random_cone
    >>> K = random_cone()
    >>> check(K)
    True
"""

from signatures import *
from cones import *


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

        >>> from sympy import expand, floor, symbols
        >>> m = symbols("m", integer=True, positive=True)
        >>> lb = lowerbound2(L(m)).replace(floor, lambda x: x)
        >>> expand(lb)
        m + 2
    """
    from signatures import f
    return 2 + f(1 + K.dim) - K.rank


def lowerbound3(K : SymmetricCone) -> int:
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
      This lower bound is always 10.

    Examples
    --------

    Along with with the other two, this one implies a lower bound of
    ``2*K.dim``. If ``n`` satisfies the other two, we can add them up
    and divide by two to get another lower bound on ``n`` in terms of
    ``K.dim``. The assumption that ``n >= 10`` takes care of cones of
    dimension five or less, and then the newly-derived inequality can
    be set greater than or equal to ``2*K.dim`` and solved to obtain a
    quadratic inequality that will be true for ``K.dim() >= 8``.
    Below we let the symbols ``d`` and ``r`` stand for ``K.dim`` and
    ``K.rank``::

        >>> from sympy import expand, floor, symbols
        >>> d,r = symbols("d,r", integer=True, positive=True)
        >>> K = SymmetricCone(0)
        >>> K.dim = d
        >>> K.rank = r
        >>> implied_bound = 2*d
        >>> g = (lowerbound1(K) + lowerbound2(K))/2 - implied_bound

    We want ``g`` to be nonnegative, but we can multiply it by ``4``
    without changing when it is nonnegative. Again we have to strip
    the symbolic `floor` ourselves because sympy doesn't know that
    `d**2 + d` is even::

        >>> g = 4*g
        >>> expand(g).replace(floor, lambda x: x)
        d**2 - 9*d + 10

    Since ``g`` is an upwards-facing parabola, it will be negative on
    an interval, and nonnegative everywhere else. We see tat for ``d >= 8``,
    g will be nonnegative::

        >>> [ g.subs({d:i}) for i in range(10) ]
        [10, 2, -4, -8, -10, -10, -8, -4, 2, 10]

    We can check the remaining two cases, ``d == 6`` and ``d == 7``,
    manually. The argument we give for Lemma 4 handles ``L(d)`` and
    all other cones as two separate cases::

        >>> from sql import all_cones_of_dim
        >>> def check(d):
        ...     lowerbound1(L(d)) >= 2*d
        ...     cs = all_cones_of_dim(d)
        ...     return all( lowerbound2(c) >= 2*d
        ...                 for c in cs
        ...                 if not c == L(d) )
        >>> check(6)
        True
        >>> check(7)
        True

    """
    return 10


def lowerbound4(K):
    r"""
    The fourth and final precondition on ``n`` that we can use
    for Theorem 4, obtained near the end of the paper by loosening
    :func:`lowerbound3`.

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
    """
    return 5
