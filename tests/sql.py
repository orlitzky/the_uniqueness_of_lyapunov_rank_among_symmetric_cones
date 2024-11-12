r"""
Functions to access the SQL database of cones and Lorentz ranks
(Lyapunov ranks arising from direct sums of Lorentz cones).

There are only two tables in the database:

    CREATE TABLE cones (
      dim INTEGER NOT NULL,
      rank INTEGER NOT NULL,
      data BLOB NOT NULL
    );
    CREATE INDEX dim_rank_idx ON cones (dim,rank);

and

    CREATE TABLE lorentz_ranks (
      dim INTEGER NOT NULL,
      rank INTEGER NOT NULL
    );
    CREATE INDEX dim_idx ON lorentz_ranks (dim);

In the "cones" table, the "data" column contains a packed, serialized
representation of the cone. In other words, the result of calling
``msgpack.packb(K.serialize())``.

Sanity check for a random dimension and rank::

    >>> import msgpack, sqlite3
    >>> from random import randint
    >>> from cones import SymmetricCone
    >>> conn = sqlite3.connect("cones.db")
    >>> cur = conn.cursor()
    >>> d = randint(0,max_cone_dim())
    >>> r = randint(d, (d**2 - d + 2//2))
    >>> stmt = "SELECT data FROM cones WHERE dim=? AND rank=?"
    >>> results = []
    >>> row = cur.execute(stmt, (d,r)).fetchone()
    >>> while row:
    ...     K = SymmetricCone.deserialize(msgpack.unpackb(row[0],
    ...                                                   use_list=False))
    ...     results.append(K.dim == d and K.rank == r)
    ...     row = cur.fetchone()
    >>> all(results)
    True
    >>> conn.close()

Ensure that the trivial cone is in the database::

    >>> conn = sqlite3.connect("cones.db")
    >>> cur = conn.cursor()
    >>> stmt = "SELECT MIN(dim) FROM cones"
    >>> cur.execute(stmt, ()).fetchone()[0] == 0
    True
    >>> conn.close()

"""
import sqlite3
import msgpack
from cones import SymmetricCone


def all_cones_of_dim(n : int) -> tuple[SymmetricCone]:
    r"""
    Return all symmetric cones having dimension ``n``.

    Setup::

    >>> from cones import *

    Base cases::

        >>> all_cones_of_dim(0)
        (L(0),)
        >>> all_cones_of_dim(1)
        (L(1),)
        >>> all_cones_of_dim(2)
        (L(1) + L(1),)

    Real examples::

        >>> sorted(all_cones_of_dim(3))
        [L(1) + L(1) + L(1), L(3)]

        >>> sorted(all_cones_of_dim(4))
        [L(1) + L(1) + L(1) + L(1), L(1) + L(3), L(4)]

        >>> sorted(all_cones_of_dim(5))
        [L(1) + L(1) + L(1) + L(1) + L(1), L(1) + L(1) + L(3), L(1) + L(4), L(5)]

    Finally at dim 6, the 3x3 real PSD cone makes an entrance::

        >>> sorted(all_cones_of_dim(6))
        [L(1) + L(1) + L(1) + L(1) + L(1) + L(1), L(1) + L(1) + L(1) + L(3), L(1) + L(1) + L(4), L(1) + L(5), L(3) + L(3), HR(3), L(6)]

    Sets of cones in different dimensions should never intersect::

        >>> from random import randint
        >>> d1 = randint(0,20)
        >>> d2 = randint(0,20)
        >>> c1 = all_cones_of_dim(d1)
        >>> c2 = all_cones_of_dim(d2)
        >>> expected = ()
        >>> if d2 == d1:
        ...     expected = c1
        >>> actual = tuple( c for c in c1 if c in c2 )
        >>> actual == expected
        True

    Check for the existence of some expected cones::

        >>> HR(3) in all_cones_of_dim(6)
        True
        >>> HC(3) in all_cones_of_dim(9)
        True
        >>> HH(3) in all_cones_of_dim(15)
        True
        >>> HO(3) in all_cones_of_dim(27)
        True

    """
    conn = sqlite3.connect("cones.db")
    cur = conn.cursor()
    stmt = "SELECT data FROM cones WHERE dim=?"
    cur.execute(stmt, (n,))
    result = tuple(
        SymmetricCone.deserialize(
          msgpack.unpackb(t[0], use_list=False)
        )
        for t in cur.fetchall()
    )
    conn.close()
    return result


def similacra(K : SymmetricCone) -> tuple[SymmetricCone]:
    r"""
    Return all similacra of the given cone.
    """
    conn = sqlite3.connect("cones.db")
    cur = conn.cursor()
    stmt = "SELECT data FROM cones WHERE dim=? AND rank=? AND data<>?"
    cur.execute(stmt, (K.dim, K.rank, msgpack.packb(K.serialize())) )
    result = tuple(
        SymmetricCone.deserialize(
          msgpack.unpackb(t[0], use_list=False)
        )
        for t in cur.fetchall()
    )
    conn.close()
    return result



def max_cone_dim():
    r"""
    Return the maximum dimension of any cone in the database.

    Examples
    --------

    This is the right answer, because I computed the database and I
    say so::

        >>> max_cone_dim()
        83

    """
    conn = sqlite3.connect("cones.db")
    cur = conn.cursor()
    stmt = "SELECT MAX(dim) FROM cones"
    result = cur.execute(stmt).fetchone()
    conn.close()
    return result[0]



def admissible_lorentz_ranks(n : int) -> tuple[int]:
    r"""
    Return the admissible Lyapunov ranks for sums of Lorentz cones
    of total dimension ``n``.

    Parameters
    ----------

    n : int
      The dimension for which you'd like to know, what Lyapunov ranks
      are possible if we consider only Lorentz cone factors?

    Returns
    -------

    tuple
      The "set" of Lyapunov ranks that can be achieved in dimension
      ``n`` using only Lorentz cone factors.

    Examples
    --------

    Low-dimensional examples::

        >>> admissible_lorentz_ranks(0)
        (0,)
        >>> admissible_lorentz_ranks(1)
        (1,)
        >>> admissible_lorentz_ranks(2)
        (2,)
        >>> sorted(admissible_lorentz_ranks(3))
        [3, 4]

    """
    conn = sqlite3.connect("cones.db")
    cur = conn.cursor()
    stmt = "SELECT rank FROM lorentz_ranks WHERE dim=?"
    cur.execute(stmt, (n,) )
    result = tuple( r[0] for r in cur.fetchall() )
    conn.close()
    return result


def admissible_ranks(n: int) -> tuple[int]:
    r"""
    Compute all admissible Lyapunov ranks for cones of total
    dimension ``n``.

    This is basically :func:`admissible_lorentz_ranks`, but with an
    optional 3x3 complex PSD factor thrown in. All symmetric cones
    share a signature with a cone of this form.

    Parameters
    ----------

    n : int
      The dimension for which you'd like to know, what Lyapunov ranks
      are possible?

    Returns
    -------

    tuple[int]
      The set of Lyapunov ranks that can be achieved in dimension ``n``

    Examples
    --------

    Low-dimensional examples agree with
    :func:`admissible_lorentz_ranks`, at least until there is room for
    that 3x3 complex PSD factor to fit::

        >>> all( admissible_ranks(k) == admissible_lorentz_ranks(k)
        ...      for k in range(9) )
        True
        >>> admissible_ranks(9) == admissible_lorentz_ranks(9)
        False

    All cones share a signature with a cone of this form::

        >>> from cones import random_cone
        >>> K = random_cone()
        >>> while K.dim > max_lorentz_rank_dim():
        ...     K = random_cone()
        >>> K.rank in admissible_ranks(K.dim)
        True

    """
    if n < 9:
        return admissible_lorentz_ranks(n)

    a = admissible_lorentz_ranks(n)
    b = tuple( 17 + r for r in admissible_lorentz_ranks(n - 9) )
    return tuple(set(a+b)) # dedupe


def max_lorentz_rank_dim():
    r"""
    Return the maximum dimension for which we know the admissible
    Lorentz ranks.

    Examples:

    This is the right answer, because I computed the database and I
    say so::

        >>> max_lorentz_rank_dim()
        250

    """
    conn = sqlite3.connect("cones.db")
    cur = conn.cursor()
    stmt = "SELECT MAX(dim) FROM lorentz_ranks"
    result = cur.execute(stmt).fetchone()
    conn.close()
    return result[0]


def insert_lorentz_ranks(n : int, ranks : list[int]):
    conn = sqlite3.connect("cones.db")
    cur = conn.cursor()
    stmt = "INSERT INTO lorentz_ranks (dim,rank) VALUES (?,?)"
    cur.executemany(stmt, ((n, r) for r in ranks) )
    conn.commit()
    conn.close()

def insert_cones(n : int, d : dict):
    conn = sqlite3.connect("cones.db")
    cur = conn.cursor()
    stmt = "INSERT INTO cones (dim,rank,data) VALUES (?,?,?)"
    cur.executemany(stmt, ((n, r, msgpack.packb(s))
                           for r in d
                           for s in d[r]) )
    conn.commit()
    conn.close()


def ranks_cones(n):
    r"""
    Base cases that should agree with :func:`dim_ranks_cones` when NOT
    using the SQL database::

        >>> ranks_cones(0)
        {0: (1,)}
        >>> ranks_cones(1)
        {1: (11,)}
        >>> ranks_cones(2)
        {2: ((11, 11),)}

    """
    conn = sqlite3.connect("cones.db")
    cur = conn.cursor()
    stmt = "SELECT rank,data FROM cones WHERE dim=?"
    result_pairs = cur.execute(stmt, (n,) ).fetchall()
    conn.close()

    # Convert the paired results to a dict
    d_n = {}
    for r,s in result_pairs:
        d_n.setdefault(r, []).append(msgpack.unpackb(s, use_list=False))

    for r in d_n:
        d_n[r] = tuple(d_n[r])

    return d_n
