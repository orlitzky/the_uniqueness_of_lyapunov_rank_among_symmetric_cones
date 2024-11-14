r"""
Exhaustive integer partition calculations to confirm results.

Test the expressions derived in Proposition 1::

    >>> from sympy import expand, floor, simplify, symbols
    >>> m,n,k = symbols("m,n,k", integer=True, positive=True)

Sympy doesn't know (for example) that n^2 + n is even, so we have to
fix the "floor" that it inserts everywhere we use integer
division-by-two on an even expression::

    >>> def fix_floor(s):
    ...     return s.replace(floor, lambda x: x)

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

Check Theorem/Conjecture 4 using our precomputed dictionary of
cones. We start with a cone ``K``, and then add ``L(n)`` factors to
it. If the resulting sum has similacra, then each similacrum should
have an ``L(n)`` factor. We do this for as many ``n`` as we can,
constrained by the fact that ``K + L(n)`` needs to have a dimension
that we have cached::

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

def lowerbound1(K):
    r"""
    The lower bound for "n" to be valid in Lemma 4, also
    appearing as the first of three lower bounds in Theorem/Conjecture
    4.

    Examples
    --------

    Check Lemma 4 using our precomputed dict of admissible Lyapunov
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
    The second lower bound for "n" in Theorem/Conjecture 4.

    Examples:

    It works with symbolic values, though we have to manually
    eliminate the `floor` that arises from python's integer division
    (m**2 + m + 2 is guaranteed to be even)::

        >>> from sympy import expand, floor, symbols
        >>> m = symbols("m", integer=True, positive=True)
        >>> lb = lowerbound2(L(m)).replace(floor, lambda x: x)
        >>> expand(lb)
        m + 2

    """
    from signatures import f
    return 2 + f(1 + K.dim) - K.rank


def lowerbound3(K):
    r"""
    The third precondition on ``n`` in Theorem (or Conjecture) 4.

    Along with with the other two, this one implies a lower bound of
    ``2*K.dim - 1``. If ``n`` satisfies the other two, we can add
    them up and divide by two to get another lower bound on ``n``. For
    ``n >= 10``, that lower bound will always be greater than
    ``2*K.dim - 1``. We can check this in sympy using ``d`` for
    ``K.dim`` and ``r`` for ``K.rank``::

        >>> from sympy import expand, floor, symbols
        >>> d,r = symbols("d,r", integer=True, positive=True)
        >>> K = SymmetricCone(0)
        >>> K.dim = d
        >>> K.rank = r
        >>> implied_bound = 2*d - 1
        >>> g = (lowerbound1(K) + lowerbound2(K))/2 - implied_bound

    We want ``g`` to be nonnegative, but we can multiply it by ``4``
    without changing when it is nonnegative. Again we have to strip
    the symbolic `floor` ourselves because sympy doesn't know that
    `d**2 + d` is even::

        >>> g = 4*g
        >>> expand(g).replace(floor, lambda x: x)
        d**2 - 9*d + 14

    Since ``g`` is an upwards-facing parabola, it will be negative on
    an interval, and nonnegative everywhere else. Here's the
    interval::

        >>> [ g.subs({d:i}) for i in range(10) ]
        [14, 6, 0, -4, -6, -6, -4, 0, 6, 14]

    """
    return 10
