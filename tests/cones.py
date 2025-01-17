r"""
Cone classes used in our test suite.
"""

# allow forward-references in type signatures
from __future__ import annotations

from itertools import chain
from math import sqrt

type SerialSum = tuple[int,...]
type SerialIrreducible = int
type SerialCone = SerialIrreducible | SerialSum

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
    # Declare attributes for the "size" n, dimension, and Lyapunov
    # rank in the superclass.
    n: int
    dim: int
    rank: int

    # All subclasses have an instance cache and an id as well.
    _instance_cache: dict
    id: int

    @staticmethod
    def irreducible_classes() -> tuple:
        return (L,HR,HC,HH,HO)

    @staticmethod
    def _deserialize_one(s : SerialIrreducible) -> SymmetricCone:
        r"""
        Deserialize one integer into an irreducible cone.

        Parameters

        s : int
          The integer to deserialize.

        Returns
        -------

        An instance of an irreducible cone (:class:`L`, :class:`HR`,
        :class:`HC`, :class:`HH`, or :class:`HO`).

        Examples
        --------

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
    def deserialize(s : SerialCone) -> SymmetricCone:
        r"""
        Deserialize either an integer or a tuple of integers into
        a symmetric cone.

        Parameters
        ----------

        s : SerialCone
          Either an int or a tuple of ints to deserialize.

        Returns
        -------

        The symmetric cone that serializes to ``s``.

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
    def name() -> str:
        r"""
        The name of this cone. This method should be overridden
        in each subclass, and exists only to make the implementation
        of :meth:`__repr__` easier.
        """
        raise NotImplementedError

    def __repr__(self) -> str:
        r"""
        The string representation of this cone.
        """
        return f"{self.name()}({self.n})"

    def __hash__(self) -> int:
        r"""
        Return a unique integer hash for this cone. The
        :meth:`serialize` method already returns a unique integer or
        tuple of integers corresponding to this cone, so we can simply
        hash that.

        WARNING: hashes in python can collide! The only guarantee
        is that equal objects must have equal hashes. If you don't
        believe this::

            >>> hash(-1) == hash(-2)
            True

        We're using hashes as the lookup keys in our instance
        caches. This should generally be safe for nonnegative integers
        (what we get when we serialize our cones), but I can't promise
        that two different tuples will never have the same hash, even
        if their components do not. So there is a small chance that we
        load the wrong :class:`DirectSum` from the instance
        cache. This would be a fun experiment to run on the cones
        database at some point: there should be no duplicate cones in
        the database, so there should be no duplicate hashes.
        """
        if self._hash is None:
            self._hash = hash(self.serialize())
        return self._hash

    def __eq__(self, other) -> bool:
        return self.serialize() == other.serialize()

    def __lt__(self, other) -> bool:
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


    @staticmethod
    def _flatten_factors(factors : tuple[SymmetricCone, ...]) -> tuple[SymmetricCone, ...]:
        r"""
        Flatten the given factors, removing all class:`DirectSum`
        wrappers.

        This is analogous to flattening a list like
        ``[a,[b,[c,d,e]]]`` down to ``[a,b,c,d,e]``.

        The implementation of method is suitable only for irreducible
        cones, and must be overridden in the class:`DirectSum` class.
        In general, this method should be suitable for use on the
        result of the class's :meth:`factors` method, which,
        for irreducible cones, will return a singleton list.

        Parameters
        ----------

        factors : tuple[SymmetricCone, ...]
          A list of symmetric cones, some of which may be direct sums.

        Returns
        -------

        A tuple of irreducible cones.

        Examples
        --------

            >>> HC(3)._flatten_factors(HC(3).factors())
            (HC(3),)

        """
        return factors


    def factors(self) -> tuple[SymmetricCone, ...]:
        r"""
        Return the factors of this symmetric cone.

        The result should be unique up to the ordering of the factors.
        This method is suitable for irreducible cones, which have only
        a single factor. It should be overridden in :class:`DirectSum`.

        To ensure that hashes are unique, we agree once and for all
        that the trivial cone has no factors.

        Returns
        -------

        A list of symmetric cones.

        TODO: we should be able to restrict the return type to
        irreducible cones.

        Examples
        --------

            >>> HC(3).factors()
            (HC(3),)

        We give the canonical answer (no factors) for the trivial
        cone::

            >>> L(0).factors()
            ()
        """
        if self.dim == 0:
            return ()
        else:
            return (self,)


    def signature(self) -> tuple[int,int]:
        r"""
        The signature of this cone.

        The signature of a cone is simply a pair consisting of its
        dimension and its Lyapunov rank.

        Returns
        -------

        A pair (2-tuple) of ints. Its first entry is the dimension of
        this cone, and its second entry is its Lypaunov rank.

        Examples
        --------

            >>> L(3).signature()
            (3, 4)

        """
        return (self.dim, self.rank)


    def serialize(self) -> SerialCone:
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


    def similacra(self) -> tuple[SymmetricCone, ...]:
        r"""
        Return all similacra of this cone.

        Returns
        -------

        A tuple of symmetric cones, not isomorphic to ``self``, but
        sharing their dimensions and Lyapunov ranks with it.

        Examples
        --------

        Examples from the paper::

            >>> HR(3).similacra()
            (L(1) + L(1) + L(4),)

        Examples from an earlier version of the paper where we computed
        similacra for multiple copies of ``HC(3)`` explicitly::

            >>> K2 = DirectSum([L(7),L(3), RN(8)])
            >>> K2 in DirectSum([HC(3)]*2).similacra()
            True
            >>> K3 = DirectSum([L(8),L(4), RN(15)])
            >>> K3 in DirectSum([HC(3)]*3).similacra()
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
    one-dimensional Lorentz cones.

    Parameters
    ----------

    n : int
      How many copies of ``L(1)`` you want.

    Returns
    -------

    The direct sum of ``n`` copies of ``L(1)``.

    Examples
    --------

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
    def name() -> str:
        return "L"

    @staticmethod
    def in_dim(d : int) -> L:
        r"""
        Return the Lorentz cone of dimension ``d``.

        This always works, because there's a Lorentz cone in every dimension.

        Examples
        --------

            >>> L.in_dim(21)
            L(21)

        """
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

# We have to fix this case because the general formula we use in
# __new__ fails at n=0.
L(0).rank = 0


class HR(SymmetricCone):
    r"""
    The real symmetric PSD cone of order ``n``.

    Examples
    --------

        >>> HR(4)
        HR(4)
        >>> HR(4).signature()
        (10, 16)

    When ``n`` is small, you will get the Lorentz cone isomorphic to
    ``HR(n)`` instead::

        >>> HR(0)
        L(0)
        >>> HR(1)
        L(1)
        >>> HR(2)
        L(3)

    """
    id = 2

    @staticmethod
    def name() -> str:
        return "HR"

    @staticmethod
    def in_dim(d : int) -> L | HR | None:
        r"""
        Return the real symmetric PSD cone of dimension ``d``, or
        ``None`` if there isn't one.

        Examples
        --------

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
    The complex hermitian PSD cone of order ``n``.

        >>> HC(4)
        HC(4)
        >>> HC(4).signature()
        (16, 31)

    When ``n`` is small, you will get the Lorentz cone isomorphic to
    ``HC(n)`` instead::

        >>> HC(0)
        L(0)
        >>> HC(1)
        L(1)
        >>> HC(2)
        L(4)

    Examples
    --------

    The direct sum of ``m >= 2`` copies of the 3-by-3 cone has
    similacra. This used to be a Lemma in the paper, but it has
    been superseded.

        >>> K2 = DirectSum([
        ...   L(7),
        ...   L(3),
        ...   RN(8)
        ... ])
        >>> K2.signature()
        (18, 34)
        >>> DirectSum(2*[HC(3)]).signature()
        (18, 34)

        >>> K3 = DirectSum([
        ...   L(8),
        ...   L(4),
        ...   RN(15)
        ... ])
        >>> K3.signature()
        (27, 51)
        >>> DirectSum(3*[HC(3)]).signature()
        (27, 51)

    Do a random check so make sure this works for many copies::

        >>> from random import randint
        >>> from sql import admissible_lorentz_ranks, max_lorentz_rank_dim
        >>> m = randint(2, 12)
        >>> K = DirectSum(m*[HC(3)])
        >>> mlrd = max_lorentz_rank_dim()
        >>> skip = (not mlrd) or (K.dim > mlrd)
        >>> skip or K.rank in admissible_lorentz_ranks(K.dim)
        True

    """
    id = 3

    @staticmethod
    def name() -> str:
        return "HC"

    @staticmethod
    def in_dim(d : int) -> L | HC | None:
        r"""
        Return the complex hermitian PSD cone of dimension ``d``,
        or ``None`` if there isn't one.

        Examples
        --------

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
    The quaternion hermitian PSD cone of order ``n``.

    Examples
    --------

        >>> HH(4)
        HH(4)
        >>> HH(4).signature()
        (28, 64)

    When ``n`` is small, you will get the Lorentz cone isomorphic to
    ``HH(n)`` instead::

        >>> HH(0)
        L(0)
        >>> HH(1)
        L(1)
        >>> HH(2)
        L(6)

    """
    id = 4

    @staticmethod
    def name() -> str:
        return "HH"

    @staticmethod
    def in_dim(d : int) -> L | HH | None:
        r"""
        Return the quaternion hermitian PSD cone of dimension ``d``, or
        ``None`` if there isn't one.

        Examples
        --------

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
    The cone of squares in the Albert algebra of order ``n``.

    Examples
    --------

    When ``n`` is small, you will get the Lorentz cone isomorphic to
    ``HO(n)`` instead::

        >>> HO(0)
        L(0)
        >>> HO(1)
        L(1)
        >>> HO(2)
        L(10)

    These cones are only symmetric for ``n <= 3``::

        >>> HO(7)
        Traceback (most recent call last):
        ...
        ValueError: invalid size n=7

    So there's only one interesting case remaining::

        >>> HO(3)
        HO(3)
        >>> HO(3).signature()
        (27, 79)

    """
    id = 5

    @staticmethod
    def name() -> str:
        return "HO"

    @staticmethod
    def in_dim(d : int) -> L | HO | None:
        r"""
        Return the octonion hermitian "PSD" cone of dimension ``d``,
        or ``None`` if there isn't one.

        Examples
        --------

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

    Examples
    --------

        >>> K = DirectSum([])
        >>> K.dim
        0
        >>> K.rank
        0

    """
    _factors: tuple[SymmetricCone, ...]

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
                f for f in cls._flatten_factors(factors)
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
    def _flatten_factors(factors : tuple[SymmetricCone, ...]) -> tuple[SymmetricCone, ...]:
        r"""
        Recursively flatten the factors of the given factors.

        The base case where we do nothing for an irreducible factor
        is implemented as :meth:`SymmetricCone._flatten_factors`.

        Parameters
        ----------

        factors : tuple[SymmetricCone, ...]
          A tuple of (not necessarily irreducible) symmetric cones.

        Returns
        -------

        A new tuple of factors, but this time with each factor
        irreducible.

        TODO: we should be able to tighten the return type!

        Examples
        --------

            >>> K = DirectSum([DirectSum([L(3)]*2), HC(3)])
            >>> K._flatten_factors(K.factors())
            (L(3), L(3), HC(3))

        """
        return sum( (f._flatten_factors(f.factors()) for f in factors),
                    () )


    def factors(self) -> tuple[SymmetricCone, ...]:
        return self._factors


    def serialize(self) -> SerialSum:
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
    Generate a random irreducible symmetric cone; that is, a
    :class:`SymmetricCone` that is not a :class:`DirectSum`.

    Returns
    -------

    An irreducible symmetric cone.

    TODO: we should be able to tighten the return type!

    Examples
    --------

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

    Returns
    -------

    A symmetric cone.

    Examples
    --------

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

    factors: list[SymmetricCone]
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

