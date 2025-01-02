r"""
Cone classes used in our test suite.
"""

from itertools import chain
from math import sqrt

# Lyapunov rank and dimension calculations
class SymmetricCone:
    r"""
    A symmetric cone of "size" ``n``.

    Equality means isomorphism in our case, and each cone should be
    represented by a unique object::

        >>> K = random_cone()
        >>> J = random_cone()
        >>> (K == J) == (J == K)
        True
        >>> (K is J) == (K == J)
        True
        >>> K == DirectSum([K]) and J == DirectSum([J])
        True

    More isomorphic examples::

        >>> RN(2) == L(2)
        True
        >>> RN(2) == DirectSum([L(1)]*2)
        True
        >>> RN(1) == DirectSum([L(1)])
        True
        >>> RN(1) == L(1)
        True
        >>> RN(1) == HR(1)
        True
        >>> RN(1) == HC(1)
        True
        >>> RN(1) == HH(1)
        True
        >>> RN(1) == HO(1)
        True
        >>> L(3) == HR(2)
        True
        >>> L(4) == HC(2)
        True
        >>> L(6) == HH(2)
        True
        >>> L(10) == HO(2)
        True
        >>> RN(4) == DirectSum([L(2)]*2)
        True

    And again. This time the isomorphism requires reordering the
    factors::

        >>> K1 = RN(3)
        >>> K2 = HR(2)
        >>> K3 = L(0)
        >>> K4 = HH(3)
        >>> K5 = L(1)
        >>> J1 = L(2)
        >>> J2 = HH(3)
        >>> J3 = L(2)
        >>> J4 = L(3)
        >>> DirectSum([K1,K2,K3,K4,K5]) == DirectSum([J1,J2,J3,J4])
        True

    """
    @staticmethod
    def irreducible_classes():
        return (L,HR,HC,HH,HO)

    @staticmethod
    def _deserialize_one(s):
        r"""
        Deserialize one integer into an irreducible cone.

            >>> SymmetricCone._deserialize_one(32)
            HR(3)
            >>> SymmetricCone._deserialize_one(29)
            Traceback (most recent call last):
            ...
            ValueError: unable to deserialize 29

        We can deserialize the trivial cone::

            >>> SymmetricCone._deserialize_one(1)
            L(0)

        """

        # Heads up: serializing L(0) produces 01 = 1.
        # The math for i and n below should produce
        # i=1 and n=0 for "1".
        i = s % 10
        n = (s - i) // 10
        for c in SymmetricCone.irreducible_classes():
            if c.id == i:
                return c(n)
        raise ValueError(f"unable to deserialize {s}")

    @staticmethod
    def deserialize(s):
        r"""
        Deserialize either an integer or a tuple of integers into
        a symmetric cone.

        Examples
        --------

        Typical examples::

            >>> SymmetricCone.deserialize(11)
            L(1)
            >>> SymmetricCone.deserialize((31, 33, 34))
            L(3) + HC(3) + HH(3)

        Invalid input (not from a serialized cone)::

            >>> SymmetricCone.deserialize((11, 100))
            Traceback (most recent call last):
            ...
            ValueError: unable to deserialize 100

        The trivial cone works as expected::

            >>> SymmetricCone.deserialize(L(0).serialize()) == L(0)
            True

        A random cone works as expected::

            >>> K = random_cone()
            >>> SymmetricCone.deserialize(K.serialize()) == K
            True

        """
        if isinstance(s, tuple):
            # The DirectSum instance cache key is actually just
            # the hash of its serialized factor tuple, i.e. what
            # we were given.
            h = hash(s)
            if h not in DirectSum._instance_cache:
                # Convert to tuple explicitly if we're going to skip
                # canonicalization.
                c = DirectSum(
                    tuple ( map(SymmetricCone._deserialize_one, s) ),
                    False
                )
                DirectSum._instance_cache[h] = c
            return DirectSum._instance_cache[h]
        else:
            return SymmetricCone._deserialize_one(s)


    def __init__(self, n):
        # cached values
        self._serial = None
        self._hash = None

        # no error checking, it makes deserialization too slow
        self.n = n

    @staticmethod
    def name():
        r"""
        The name of this cones. This method should be overridden
        in each subclass, and exists only to make the implementation
        of :meth:`__repr__` easier.
        """
        raise NotImplementedError

    def __repr__(self):
        r"""
        The string representation of this cone.
        """
        return f"{self.name()}({self.n})"

    def __hash__(self):
        r"""
        Return a unique integer hash for this cone. The
        :meth:`serialize` method already returns a unique integer or
        tuple of integers corresponding to this cone, so we can simply
        hash that.
        """
        if self._hash is None:
            self._hash = hash(self.serialize())
        return self._hash

    def __eq__(self, other):
        return self.serialize() == other.serialize()

    def __lt__(self, other):
        r"""
        Implement an ordering for symmetric cones.

        We already have equality via :meth:`serialize`, and we can
        derive "less than" from that too. Basically we treat
        irreducible cones as singleton direct sums, serialize
        everything to tuples of integers, and then use the
        lexicographic "less than" for tuples.

        Examples
        --------

        A simple example where "less than" is consistent::

            >>> L(3) < L(3)
            False
            >>> L(3) < L(4)
            True
            >>> L(4) < L(3)
            False

        Larger cones of the same type should be greater than smaller
        cones. One nice aspect of the integer sort is that it handles
        this correctly whereas an alphabetical sort of the dimension
        would not::

            >>> L(1) < L(10)
            True
            >>> L(1) < L(2)
            True
            >>> L(10) < L(2)
            False

        Another pleasing side effect of our :meth:`serialize` method
        is that it sorts the real, complex, quaternion, and octonion
        PSD cones of the same size all in that order::

            >>> from random import randint
            >>> n = randint(2,10)
            >>> HR(n) < HC(n) < HH(n)
            True
            >>> n > 3 or  HH(n) < HO(n)
            True

        This implementation of "less than" leads to a total order::

            >>> K = random_irreducible_cone()
            >>> J = random_irreducible_cone()
            >>> (K < J) or (J < K) or (K == J)
            True
            >>> not ((K < J) and (J < K))
            True
            >>> (K != J) or not (K < J or J < K)
            True
            >>> H = random_irreducible_cone()
            >>> (not (J < K and H < J)) or (H < K)
            True
            >>> (not (K < J and J < H)) or (K < H)
            True

        """
        s = self.serialize()
        o = other.serialize()
        if not isinstance(s, tuple):
            s = (s,)
        if not isinstance(o, tuple):
            o = (o,)

        return s < o

    def signature(self):
        r"""
        The signature of this cone.

        The signature of a cone is simply a pair consisting of its
        dimension and its Lyapunov rank.
        """
        return (self.dim, self.rank)

    def serialize(self):
        r"""
        Serialize this cone to an integer (irreducible cones) or
        tuple of integers (direct sums).

        No two nonequal cones should serialize to the same value. The
        default is appropriate only for irreducible cones.

        Examples
        --------

        The one unusual case is the trivial cone, where ``01`` becomes
        ``1``, but this should not cause any problems so long as we are
        expecting it::

            >>> L(0).serialize()
            1

        """
        if self._serial is None:
            self._serial = 10*self.n + self.id
        return self._serial

    def similacra(self) -> tuple:
        r"""
        Return all similacra of this cone.

        Examples from the paper::

            >>> HR(3).similacra()
            (L(1) + L(1) + L(4),)

        Lemma 1::

            >>> K2 = DirectSum([L(7),L(3), RN(8)])
            >>> K2 in DirectSum([HC(3)]*2).similacra()
            True
            >>> K3 = DirectSum([L(8),L(4), RN(15)])
            >>> K3 in DirectSum([HC(3)]*3).similacra()
            True

        Proposition 5::

            >>> K = DirectSum([L(5), L(5), L(4), RN(2)])
            >>> K in HC(4).similacra()
            True

        Proposition 6::

            >>> K = DirectSum([L(8), RN(7)])
            >>> K in HH(3).similacra()
            True

            >>> K = DirectSum([L(10), RN(18)])
            >>> K in HH(4).similacra()
            True

            >>> K = DirectSum([L(12), RN(33)])
            >>> K in HH(5).similacra()
            True

        Proposition 7::

            >>> K = DirectSum([L(11), L(5), L(3), RN(8)])
            >>> K in HO(3).similacra()
            True

        A cone is never its own similacrum::

            >>> import sql
            >>> K = random_cone()
            >>> while K.dim > sql.max_cone_dim():
            ...     K = random_cone()
            >>> K in K.similacra()
            False

        """
        from sql import similacra
        return similacra(self)


def RN(n : int) -> SymmetricCone:
    r"""
    A convenient wrapper around an ``n``-fold direct sum of
    one-dimensional Lorentz cones::

        >>> RN(3)
        L(1) + L(1) + L(1)
        >>> RN(0)
        L(0)

    """
    return DirectSum((L(1),)*n)


class L(SymmetricCone):
    r"""
    The Lorentz cone in dimension ``n``.

    Examples::

        >>> L(3)
        L(3)

    """
    id = 1

    @staticmethod
    def name():
        return "L"

    @staticmethod
    def in_dim(d):
        return L(d)

    # Instance cache. After the class definition, we'll prepopulate
    # it with L(0) and fix its Lyapunov rank so we can avoid special
    # cases in the constructor.
    _instance_cache = {}
    def __new__(cls, n):
        r"""
        Examples::

            >>> L(3).dim
            3
            >>> L(3).rank
            4

            >>> L(3) == L(3)
            True
            >>> L(3) is L(3)
            True

        And the edge case that the usual formula fails to account for::

            >>> L(0).rank
            0

        """
        # there's a special case for n=2 inserted after the DirectSum
        # class is defined
        if not n in cls._instance_cache:
            cls._instance_cache[n] = super().__new__(cls)
            cls._instance_cache[n].dim = n
            cls._instance_cache[n].rank = (n**2 - n + 2)//2
        return cls._instance_cache[n]

# fix this, the general formula is wrong.
L(0).rank = 0


class HR(SymmetricCone):
    r"""
    The real PSD cone of order ``n``.

    Examples
    --------

    Test the inverse dimension formula::

        >>> HR.in_dim(0)
        L(0)
        >>> HR.in_dim(1)
        L(1)
        >>> HR.in_dim(2)
        >>> HR.in_dim(3)
        L(3)
        >>> HR.in_dim(4)
        >>> HR.in_dim(5)
        >>> HR.in_dim(6)
        HR(3)
        >>> HR.in_dim(7)
        >>> HR.in_dim(8)
        >>> HR.in_dim(9)
        >>> HR.in_dim(10)
        HR(4)

    """
    id = 2

    @staticmethod
    def name():
        return "HR"

    @staticmethod
    def in_dim(d):
        n = (sqrt(8*d + 1) - 1)/2
        if n.is_integer():
            return HR(int(n))
        else:
            return None

    # Instance cache. All zero- or one-dimensional cones are L(0) or
    # L(1). The 2x2 cone is canonically isomorphic to the 3d Lorentz
    # cone.
    _instance_cache = { 0: L(0), 1: L(1), 2: L(3) }
    def __new__(cls, n):
        if not n in cls._instance_cache:
            cls._instance_cache[n] = super().__new__(cls)
            cls._instance_cache[n].dim = (n**2 + n) // 2
            cls._instance_cache[n].rank = n**2
        return cls._instance_cache[n]


class HC(SymmetricCone):
    r"""
    The direct sum of ``m >= 2`` copies of the 3-by-3 cone has
    similacra (Lemma 1)::

        >>> K2 = DirectSum([
        ...   L(7),
        ...   L(3),
        ...   RN(8)
        ... ])
        >>> K3 = DirectSum([
        ...   L(8),
        ...   L(4),
        ...   RN(15)
        ... ])
        >>> K2.signature()
        (18, 34)
        >>> H = HC(3)
        >>> DirectSum([H,H]).signature()
        (18, 34)
        >>> K3.signature()
        (27, 51)
        >>> DirectSum([H,H,H]).signature()
        (27, 51)

    Do a random check so make sure this works for many copies::

        >>> from random import randint
        >>> from sql import admissible_lorentz_ranks
        >>> m = randint(2, 12)
        >>> K = DirectSum(m*[H])
        >>> K.rank in admissible_lorentz_ranks(K.dim)
        True

    Test the inverse dimension formula::

        >>> HC.in_dim(0)
        L(0)
        >>> HC.in_dim(1)
        L(1)
        >>> HC.in_dim(2)
        >>> HC.in_dim(3)
        >>> HC.in_dim(4)
        L(4)
        >>> HC.in_dim(5)
        >>> HC.in_dim(6)
        >>> HC.in_dim(7)
        >>> HC.in_dim(8)
        >>> HC.in_dim(9)
        HC(3)
        >>> HC.in_dim(10)
        >>> HC.in_dim(11)
        >>> HC.in_dim(12)
        >>> HC.in_dim(13)
        >>> HC.in_dim(14)
        >>> HC.in_dim(15)
        >>> HC.in_dim(16)
        HC(4)

    """
    id = 3

    @staticmethod
    def name():
        return "HC"

    @staticmethod
    def in_dim(d):
        n = sqrt(d)
        if n.is_integer():
            return HC(int(n))
        else:
            return None

    # Instance cache. All zero- or one-dimensional cones are L(0) or
    # L(1). The 2x2 cone is canonically isomorphic to the 4d Lorentz
    # cone.
    _instance_cache = { 0: L(0), 1: L(1), 2: L(4) }
    def __new__(cls, n):
        if not n in cls._instance_cache:
            cls._instance_cache[n] = super().__new__(cls)
            cls._instance_cache[n].dim = n**2
            cls._instance_cache[n].rank = 2*n**2 - 1
        return cls._instance_cache[n]


class HH(SymmetricCone):
    r"""

    We know the formula for each similacrum from Proposition 5. First
    we check ``n == 3``, ``n == 4``, and ``n == 5`` individually::


    Now we check the remaining ``n >= 6`` using the generic formula::

        >>> def check(n):
        ...     K1 = HC(n+1)
        ...     K2 = L(n+1)
        ...     K3 = L(n+1)
        ...     K4 = RN(n**2 - 5*n - 3)
        ...     K = DirectSum([K1,K2,K3,K4])
        ...     return (HH(n).signature() == K.signature())
        >>>
        >>> all( check(n) for n in range(6,100) )
        True

    Test the inverse dimension formula::

        >>> HH.in_dim(0)
        L(0)
        >>> HH.in_dim(1)
        L(1)
        >>> HH.in_dim(2)
        >>> HH.in_dim(3)
        >>> HH.in_dim(4)
        >>> HH.in_dim(5)
        >>> HH.in_dim(6)
        L(6)
        >>> HH.in_dim(7)
        >>> HH.in_dim(8)
        >>> HH.in_dim(9)
        >>> HH.in_dim(10)
        >>> HH.in_dim(11)
        >>> HH.in_dim(12)
        >>> HH.in_dim(13)
        >>> HH.in_dim(14)
        >>> HH.in_dim(15)
        HH(3)

    """
    id = 4

    @staticmethod
    def name():
        return "HH"

    @staticmethod
    def in_dim(d):
        if d == 0:
            # This is the one special case where the other solution to
            # the quadratic formula is valid.
            return HH(0)

        n = (sqrt(8*d + 1) + 1)/4
        if n.is_integer():
            return HH(int(n))
        else:
            return None

    # Instance cache. All zero- or one-dimensional cones are L(0) or
    # L(1). The 2x2 cone is canonically isomorphic to the 6d Lorentz
    # cone.
    _instance_cache = { 0: L(0), 1: L(1), 2: L(6) }
    def __new__(cls, n):
        if not n in cls._instance_cache:
            cls._instance_cache[n] = super().__new__(cls)
            cls._instance_cache[n].dim = 2*n**2 - n
            cls._instance_cache[n].rank = 4*n**2
        return cls._instance_cache[n]


class HO(SymmetricCone):
    r"""
    There's only one of these (up to isomorphism) and it has a
    similacrum::

        >>> K1 = L(11)
        >>> K2 = L(5)
        >>> K3 = L(3)
        >>> K4 = RN(8)
        >>> K = DirectSum([K1,K2,K3,K4])
        >>> K.signature() == HO(3).signature()
        True

    These cones are only symmetric for ``n <= 3``::

        >>> HO(7)
        Traceback (most recent call last):
        ...
        ValueError: invalid size n=7

    Test the inverse dimension formula::

        >>> HO.in_dim(0)
        L(0)
        >>> HO.in_dim(1)
        L(1)
        >>> HO.in_dim(2)
        >>> HO.in_dim(4)
        >>> HO.in_dim(3)
        >>> HO.in_dim(5)
        >>> HO.in_dim(6)
        >>> HO.in_dim(7)
        >>> HO.in_dim(8)
        >>> HO.in_dim(9)
        >>> HO.in_dim(10)
        L(10)
        >>> HO.in_dim(26)
        >>> HO.in_dim(27)
        HO(3)

    """
    id = 5

    @staticmethod
    def name():
        return "HO"

    @staticmethod
    def in_dim(d):
        r"""
        The would-be Octonion cone in dimension 52 is not symmetric
        (it corresponds to n=4)::

            >>> HO.in_dim(52) is None
            True

        """
        # why bother with the quadratic formula when we know all of them?
        if d == 0:
            return L(0)
        elif d == 1:
            return L(1)
        elif d == 10:
            return L(10)
        elif d == 27:
            return HO(3)
        else:
            return None

    # Instance cache. All zero- or one-dimensional cones are L(0) or
    # L(1). The 2x2 cone is canonically isomorphic to the 10d Lorentz
    # cone.
    _instance_cache = { 0: L(0), 1: L(1), 2: L(10) }
    def __new__(cls, n):
        r"""
        Examples::

        >>> HO(0).dim
        0
        >>> HO(1).dim
        1
        >>> HO(2).dim
        10
        >>> HO(3).dim
        27
        """
        if not n in cls._instance_cache:
            if n > 3:
                raise ValueError(f"invalid size n={n}")
            cls._instance_cache[n] = super().__new__(cls)
            cls._instance_cache[n].dim = 27
            cls._instance_cache[n].rank = 79
        return cls._instance_cache[n]


class DirectSum(SymmetricCone):
    r"""
    A direct sum of factors, provided as a generator (list,
    tuple, set, whatever).

    Examples::

        >>> K = DirectSum([])
        >>> K.dim
        0
        >>> K.rank
        0

    """

    # instance cache
    _instance_cache = {}
    def __new__(cls, factors, canonicalize : bool = True):
        r"""
        Ensure that factors are flattened, and that trivial cones
        are removed::

            >>> DirectSum([
            ...   L(2),
            ...   DirectSum([ L(0), HR(3) ]),
            ...   DirectSum([ HC(3),
            ...               DirectSum([HO(3), L(1)])
            ...   ])
            ... ])
            L(1) + L(1) + L(1) + HR(3) + HC(3) + HO(3)

            >>> DirectSum([HC(3),L(0)]) == HC(3)
            True

        """
        # We have to remove trivial factors before we do anything
        # else, otherwise the length of the factors list might be
        # wrong. For example we need [HC(3),L(0)] to return HC(3),
        # which only happens if there is one factor.
        factors = [f for f in factors if not f.dim == 0]
        lf = len(factors)
        if lf == 0:
            return L(0)
        elif lf == 1:
            return factors[0]

        # len(factors) >= 2 guaranteed
        #
        # there's no "n" for a direct sum, but the hash of our
        # serialization is eventually what the hash of this cone
        # will be, and a hash is exactly what we want for the key
        if canonicalize:
            # we'll get the wrong hash if we don't flatten first!
            factors = tuple(sorted(
                f for f in DirectSum._flatten_factors(factors)
            ))

        _serial = tuple( f.serialize() for f in factors )
        _hash = hash(_serial)
        if not _hash in cls._instance_cache:
            K = super().__new__(cls)
            # just computed these, might as well save them
            K._factors = factors
            K._serial = _serial
            K._hash = _hash
            K.dim = 0
            K.rank = 0
            for f in factors:
                K.dim += f.dim
                K.rank += f.rank
            cls._instance_cache[_hash] = K
        return cls._instance_cache[_hash]

    def __init__(self, factors, canonicalize : bool = True):
        pass

    def __repr__(self):
        return " + ".join(repr(f) for f in self._factors)

    @staticmethod
    def _flatten_factors(factors):
        r"""
        Flatten an iterable of factors by removing all ``DirectSum``
        wrappers. This is analogous to flattening a list like
        ``[a,[b,[c,d,e]]]`` down to ``[a,b,c,d,e]``.
        """
        result = []
        for f in factors:
            if isinstance(f, DirectSum):
                result += DirectSum._flatten_factors(f._factors)
            else:
                result.append(f)

        return result

    def factors(self):
        return self._factors

    def serialize(self):
        r"""
        Serialize this cone to a tuple of integers.

        We override the superclass default to return the tuple we get
        from serializing our factors.
        """
        if self._serial is None:
            self._serial = tuple( f.serialize() for f in self._factors )
        return self._serial

# The 2d Lorentz cone is canonically isomorphic
# to a pair of 1d Lorentz cones.
L._instance_cache[2] = DirectSum( (L(1),L(1)) )

def random_irreducible_cone() -> SymmetricCone:
    r"""
    Generate a random irreducible (i.e. not a :class:`DirectSum`)
    symmetric cone.

    Examples::

        >>> isinstance(random_irreducible_cone(), DirectSum)
        False

    """
    from random import choice, randint

    c = choice(SymmetricCone.irreducible_classes())
    n = randint(1,10)
    if c == HO:
        # otherwise not symmetric
        n = min(3,n)
    elif c == L and n == 2:
        # otherwise not irreducible... return our friend?
        c = HC
        n = 3
    return c(n)


def random_cone() -> SymmetricCone:
    r"""
    Produce a somewhat-random cone with ten or fewer factors.

    Examples::

        >>> K = random_cone()
        >>> K.dim >= 0
        True
        >>> K.rank >= K.dim
        True
        >>> isinstance(K, SymmetricCone)
        True

    """
    from random import choice, randint

    # Produce at most 10 factors...
    num_factors = randint(1, 10)

    factors = []

    # All of the others share signatures with Lorentz cones, so we
    # bump up the probability that a 3x3 complex PSD factor will get
    # selected here.
    if randint(1,4) == 4:
        factors.append(HC(3))
        num_factors -= 1

    classes = [RN, L, HR, HC, HH, HO]

    while (len(factors) < num_factors):
        # ... with each having at most "size" 6.
        n = randint(1, 10)
        c = choice(classes)

        if c == HO:
            # This isn't symmetric for n > 3.
            n = min(3,n)

        factors.append(c(n))

    return DirectSum(factors)
