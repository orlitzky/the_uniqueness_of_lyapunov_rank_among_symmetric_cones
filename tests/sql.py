r"""
Functions to access the SQL database of cones and Lorentz ranks
(Lyapunov ranks arising from direct sums of Lorentz cones).

There are only two tables in the database, "cones" and
"lorentz_ranks". Their structure can be inferred from
:func:`new_database`.

Sanity check for a random dimension and rank::

    >>> import sqlite3
    >>> from msgpack import unpackb
    >>> from random import randint
    >>> from cones import SymmetricCone
    >>> d = randint(0, max_cone_dim())
    >>> r = randint(d, (d**2 - d + 2//2))
    >>> conn = sqlite3.connect(LIVE_DATABASE)
    >>> stmt = "SELECT data FROM cones WHERE dim=? AND rank=?"
    >>> result = False
    >>> with conn:
    ...     result = all(
    ...       K.dim == d and K.rank == r
    ...       for row in conn.execute(stmt, (d,r)).fetchall()
    ...       if (K := SymmetricCone.deserialize(unpackb(row[0],
    ...                                          use_list=False)))
    ...     )
    >>> result
    True
    >>> conn.close()

Ensure that the trivial cone is in the database::

    >>> conn = sqlite3.connect(LIVE_DATABASE)
    >>> stmt = "SELECT MIN(dim) FROM cones"
    >>> with conn:
    ...     conn.execute(stmt, ()).fetchone()[0] == 0
    True
    >>> conn.close()
"""
import sqlite3
import msgpack
from cones import SymmetricCone

# The default "live" database.
LIVE_DATABASE = "cones.db"

# The database used for testing.
TEST_DATABASE = "tests.db"

def new_database(db : str = TEST_DATABASE):
    r"""
    Create a new database (destroying any existing databases of
    the given name) containing our two tables.

    This is destructive, so we use the test database as the default.

    Parameters
    ----------

    db : str, default=TEST_DATABASE
      The name of the SQLite database to use.

    Examples
    --------

    Create a fresh test database::

        >>> new_database()
        >>> conn = sqlite3.connect(TEST_DATABASE)
        >>> stmt = "SELECT name FROM sqlite_master WHERE type='table';"
        >>> with conn:
        ...     print(conn.execute(stmt).fetchall())
        [('cones',), ('lorentz_ranks',), ('partitions',)]
        >>> conn.close()

    """
    from os import remove
    try:
        remove(db)
    except FileNotFoundError:
        pass

    conn = sqlite3.connect(db)

    with conn:
        conn.execute("""CREATE TABLE cones (
          dim INTEGER NOT NULL,
          rank INTEGER NOT NULL,
          data BLOB NOT NULL
        );""")
        conn.execute("CREATE INDEX dim_rank_idx ON cones (dim,rank);")

        conn.execute("""CREATE TABLE lorentz_ranks (
          dim INTEGER NOT NULL,
          rank INTEGER NOT NULL
        );""")
        conn.execute("CREATE INDEX dim_idx ON lorentz_ranks (dim);")

        # one partition (msgpack'd list) per row
        conn.execute("""CREATE TABLE partitions (
          n INTEGER NOT NULL,
          partition BLOB NOT NULL
        );""")
        conn.execute("CREATE INDEX n_idx ON partitions (n);")

    conn.close()


def partitions_of(n : int, db : str = LIVE_DATABASE) -> list[list[int]]:
    r"""
    Return all partitions of the integer ``n``.

    Since this does not modify the database, we use the live database
    by default.

    Parameters
    ----------

    n : int
      The integer you want to partition.

    db : str, default=LIVE_DATABASE
      The name of the SQLite database to use.

    Examples
    --------

        >>> new_database(db=TEST_DATABASE)
        >>> insert_partitions(0, [[0]])
        >>> insert_partitions(1, [[1]])
        >>> insert_partitions(2, [[1,1], [2]])
        >>> partitions_of(0, db=TEST_DATABASE)
        [[0]]
        >>> partitions_of(1, db=TEST_DATABASE)
        [[1]]
        >>> partitions_of(2, db=TEST_DATABASE)
        [[1, 1], [2]]

    """
    conn = sqlite3.connect(db)
    stmt = "SELECT partition FROM partitions WHERE n=?"
    result = []
    with conn:
        result = [
            msgpack.unpackb(t[0])
            for t in conn.execute(stmt, (n,)).fetchall()
        ]
    conn.close()
    return result


def insert_partitions(n : int, ps : list[list[int]], db : str = TEST_DATABASE):
    r"""
    Insert one or more partitions of ``n`` into the database.

    This is destructive, so we use the test database as the default.

    Parameters
    ----------

    n : int
      The integer that has been partitioned.

    ps : list[list[int]]
      A list of partitions, where each partition consists of a list of
      integers that sum to ``n``.

    db : str, default=TEST_DATABASE
      The name of the SQLite database to use.

    Returns
    -------

    Nothing.

    Examples
    --------

    A simple example::

        >>> new_database(db=TEST_DATABASE)
        >>> insert_partitions(2, [[3], [4,5]])
        >>> partitions_of(2, db=TEST_DATABASE)
        [[3], [4, 5]]

    """
    conn = sqlite3.connect(db)
    stmt = "INSERT INTO partitions (n,partition) VALUES (?,?)"
    with conn:
        conn.executemany(stmt, ((n, msgpack.packb(p)) for p in ps) )
    conn.close()


def all_cones_of_dim(n : int, db : str = LIVE_DATABASE) -> tuple[SymmetricCone]:
    r"""
    Return all symmetric cones having dimension ``n``.

    Since this does not modify the database, we use the live database
    by default.

    Parameters
    ----------

    n : int
      The dimension of the cones you want.

    db : str, default=LIVE_DATABASE
      The name of the SQLite database to use.

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
    conn = sqlite3.connect(db)
    stmt = "SELECT data FROM cones WHERE dim=?"
    result = ()
    with conn:
        result = tuple(
            SymmetricCone.deserialize(
              msgpack.unpackb(t[0], use_list=False)
            )
            for t in conn.execute(stmt, (n,)).fetchall()
        )
    conn.close()
    return result


def similacra(K : SymmetricCone, db : str = LIVE_DATABASE) -> tuple[SymmetricCone]:
    r"""
    Return all similacra of the given cone.

    Since this does not modify the database, we use the live database
    by default.

    Parameters
    ----------

    K : SymmetricCone
      The cone whose similacra you want.

    db : str, default=LIVE_DATABASE
      The name of the SQLite database to use.

    """
    conn = sqlite3.connect(db)
    stmt = "SELECT data FROM cones WHERE dim=? AND rank=? AND data<>?"
    args = ( K.dim, K.rank, msgpack.packb(K.serialize()) )
    result = ()
    with conn:
        result = tuple(
            SymmetricCone.deserialize(
              msgpack.unpackb(t[0], use_list=False)
            )
            for t in conn.execute(stmt, args).fetchall()
        )
    conn.close()
    return result



def have_cone_dim(d : int, db : str = LIVE_DATABASE) -> bool:
    r"""
    Return ``True`` if we have cached the cones of the given
    dimension, and ``False`` otherwise.

    Since this does not modify the database, we use the live database
    by default. This is _almost_ redundant, since we should have all
    cones of dimension less than func:`max_cone_dim`, but it's nicer
    to check for exactly what we need.

    Parameters
    ----------

    d : int
      The dimension to check for in the database.

    db : str, default=LIVE_DATABASE
      The name of the SQLite database to use.

    Examples
    --------

    An outrageously large example::

        >>> have_cone_dim(8675309)
        False

    In a new database, we won't have anything...::

       >>> new_database(db=TEST_DATABASE)
       >>> all( not have_cone_dim(n, db=TEST_DATABASE)
       ...      for n in range(25) )
       True

    ...until we add it::

       >>> import compute
       >>> _ = compute.dim_ranks_cones(5, True)
       >>> have_cone_dim(5, db=TEST_DATABASE)
       True

    """
    conn = sqlite3.connect(db)
    stmt = "SELECT dim FROM cones where dim=?"
    result = False
    with conn:
        if conn.execute(stmt, (d,)).fetchone():
            result = True
    conn.close()
    return result


def max_cone_dim(db : str = LIVE_DATABASE) -> int:
    r"""
    Return the maximum dimension of any cone in the database.

    Since this does not modify the database, we use the live database
    by default.

    Parameters
    ----------

    db : str, default=LIVE_DATABASE
      The name of the SQLite database to use.

    Examples
    --------

    The right answer depends on how long you're willing to wait (and
    whether or not the data exist)::

        >>> mcd = max_cone_dim()
        >>> mcd is None or isinstance(mcd, int)
        True

    In a new database, there won't be a maximum::

       >>> new_database(db=TEST_DATABASE)
       >>> print(max_cone_dim(db=TEST_DATABASE))
       None

    """
    conn = sqlite3.connect(db)
    stmt = "SELECT MAX(dim) FROM cones"
    result = 0
    with conn:
        result = conn.execute(stmt).fetchone()[0]
    conn.close()
    return result



def admissible_lorentz_ranks(n : int, db : str = LIVE_DATABASE) -> tuple[int]:
    r"""
    Return the admissible Lyapunov ranks for sums of Lorentz cones
    of total dimension ``n``.

    Since this does not modify the database, we use the live database
    by default.

    Parameters
    ----------

    n : int
      The dimension for which you'd like to know: what Lyapunov ranks
      are possible if we consider only Lorentz cone factors?

    db : str, default=LIVE_DATABASE
      The name of the SQLite database to use.

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
    conn = sqlite3.connect(db)
    stmt = "SELECT rank FROM lorentz_ranks WHERE dim=?"
    result = ()
    with conn:
        result = tuple( r[0]
                        for r in conn.execute(stmt, (n,)).fetchall() )
    conn.close()
    return result


def admissible_ranks(n: int, db : str = LIVE_DATABASE) -> tuple[int]:
    r"""
    Compute all admissible Lyapunov ranks for cones of total
    dimension ``n``.

    This is basically :func:`admissible_lorentz_ranks`, but with an
    optional 3x3 complex PSD factor thrown in. All symmetric cones
    share a signature with a cone of this form.

    Since this does not modify the database, we use the live database
    by default.

    Parameters
    ----------

    n : int
      The dimension for which you'd like to know: what Lyapunov ranks
      are possible?

    db : str, default=LIVE_DATABASE
      The name of the SQLite database to use.

    Returns
    -------

    tuple[int]
      The set of Lyapunov ranks that can be achieved in dimension ``n``

    db : str, default=TEST_DATABASE
      The name of the SQLite database to use.

    Examples
    --------

    Low-dimensional examples agree with
    :func:`admissible_lorentz_ranks`, at least until there is room for
    that 3x3 complex PSD factor to fit::

        >>> all( admissible_ranks(k)
        ...      ==
        ...      admissible_lorentz_ranks(k)
        ...      for k in range(9) )
        True
        >>> ( admissible_ranks(9)
        ...   ==
        ...   admissible_lorentz_ranks(9) )
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
        return admissible_lorentz_ranks(n, db)

    a = admissible_lorentz_ranks(n, db)
    b = tuple( 17 + r for r in admissible_lorentz_ranks(n - 9, db) )
    return tuple(set(a+b)) # dedupe



def have_lorentz_rank_dim(d : int, db : str = LIVE_DATABASE) -> bool:
    r"""
    Return ``True`` if we know the admissible Lorentz ranks for
    the given dimension, and ``False`` otherwise.

    Since this does not modify the database, we use the live database
    by default. This is _almost_ redundant, since we should have all
    ranks for dimensions less than func:`max_lorentz_rank_dim`, but
    it's nicer to check for exactly what we need.

    Parameters
    ----------

    d : int
      The dimension to check for in the database.

    db : str, default=LIVE_DATABASE
      The name of the SQLite database to use.

    Examples
    --------

    An outrageously large example::

        >>> have_lorentz_rank_dim(8675309)
        False

    In a new database, we won't have anything...::

       >>> new_database(db=TEST_DATABASE)
       >>> all( not have_lorentz_rank_dim(n, db=TEST_DATABASE)
       ...      for n in range(25) )
       True

    ...until we add it::

       >>> import compute
       >>> _ = compute.admissible_lorentz_ranks(10, True)
       >>> have_lorentz_rank_dim(10, db=TEST_DATABASE)
       True
    """
    conn = sqlite3.connect(db)
    stmt = "SELECT dim FROM lorentz_ranks where dim=?"
    result = False
    with conn:
        if conn.execute(stmt, (d,)).fetchone():
            result = True
    conn.close()
    return result


def max_lorentz_rank_dim(db : str = LIVE_DATABASE) -> int:
    r"""
    Return the maximum dimension for which we know the admissible
    Lorentz ranks.

    Since this does not modify the database, we use the live database
    by default.

    Parameters
    ----------

    db : str, default=LIVE_DATABASE
      The name of the SQLite database to use.

    Examples
    --------

    The right answer depends on how long you're willing to wait (and
    whether or not the data exist)::

        >>> mlrd = max_lorentz_rank_dim()
        >>> mlrd is None or isinstance(mlrd, int)
        True

    In a new database, there won't be a maximum::

       >>> new_database(db=TEST_DATABASE)
       >>> print(max_lorentz_rank_dim(db=TEST_DATABASE))
       None

    """
    conn = sqlite3.connect(db)
    stmt = "SELECT MAX(dim) FROM lorentz_ranks"
    result = 0
    with conn:
        result = conn.execute(stmt).fetchone()[0]
    conn.close()
    return result


def insert_lorentz_ranks(n : int, ranks : list[int], db : str = TEST_DATABASE):
    r"""
    Insert one dimension's worth of admissible Lorentz ranks into
    the database.

    This is destructive, so we use the test database as the default.

    Parameters
    ----------

    n : int
      The dimension to which your list of Lyapunov ranks corresponds.

    ranks : list[int]
      A list of Lyapunov ranks that are achievable using only Lorentz
      cones in dimension ``n``.

    db : str, default=TEST_DATABASE
      The name of the SQLite database to use.

    Returns
    -------

    Nothing.

    Examples
    --------

    A simple example::

        >>> new_database(db=TEST_DATABASE)
        >>> insert_lorentz_ranks(8, [6,7,5,3,0,9])
        >>> stmt = "SELECT dim,rank FROM lorentz_ranks"
        >>> conn = sqlite3.connect(TEST_DATABASE)
        >>> with conn:
        ...     print(conn.execute(stmt).fetchall())
        [(8, 6), (8, 7), (8, 5), (8, 3), (8, 0), (8, 9)]
        >>> conn.close()

    """
    conn = sqlite3.connect(db)
    stmt = "INSERT INTO lorentz_ranks (dim,rank) VALUES (?,?)"
    with conn:
        conn.executemany(stmt, ((n, r) for r in ranks) )
    conn.close()


def insert_cones(n : int, d : dict, db : str = TEST_DATABASE):
    r"""
    Insert cones into the database from a rank => cones dictionary.

    This is used in the implementation of :func:`compute.dim_ranks_cones`
    to save newly-computed cones.

    This is destructive, so we use the test database as the default.

    Parameters
    ----------

    n : int
      The dimension of your cones.

    d : dict
      A dictionary whose keys are Lyapunov ranks and whose values
      are tuples consisting of all serialized cones of dimension
      ``n`` having those Lyapunov ranks.

    db : str, default=TEST_DATABASE
      The name of the SQLite database to use.

    Returns
    -------

    Nothing.

    Examples
    --------

    Insert the two cones of dimension three, and then check that we
    can pull them back out with :func:`dim_ranks_cones`::

        >>> from cones import L, RN
        >>> new_database(db=TEST_DATABASE)
        >>> n = 3
        >>> d = { 3: [RN(3).serialize()], 4: [L(3).serialize()] }
        >>> insert_cones(n, d, db=TEST_DATABASE)
        >>> dim_ranks_cones(3, db=TEST_DATABASE)
        {3: ((11, 11, 11),), 4: (31,)}

    """
    conn = sqlite3.connect(db)
    stmt = "INSERT INTO cones (dim,rank,data) VALUES (?,?,?)"
    with conn:
        conn.executemany(stmt, ((n, r, msgpack.packb(s))
                                for r in d
                                for s in d[r]) )
    conn.close()


def dim_ranks_cones(n : int, db : str = LIVE_DATABASE):
    r"""
    Return all cones of dimension ``n`` as a rank => cones map.

    Since this does not modify the database, we use the live database
    by default.

    Parameters
    ----------

    n : int
      The dimension of the cones you want.

    db : str, default=LIVE_DATABASE
      The name of the SQLite database to use.

    Returns
    -------

    A dictionary whose keys are Lyapunov ranks and whose values
    are tuples consisting of all serialized cones of dimension
    ``n`` having those Lyapunov ranks.

    Examples
    --------

    Base cases that should agree with :func:`compute.dim_ranks_cones` when NOT
    using the SQL database::

        >>> dim_ranks_cones(0)
        {0: (1,)}
        >>> dim_ranks_cones(1)
        {1: (11,)}
        >>> dim_ranks_cones(2)
        {2: ((11, 11),)}

    """
    conn = sqlite3.connect(db)
    stmt = "SELECT rank,data FROM cones WHERE dim=?"
    result_pairs = []
    with conn:
        result_pairs = conn.execute(stmt, (n,) ).fetchall()
    conn.close()

    # Convert the paired results to a dict
    d_n = {}
    for r,s in result_pairs:
        d_n.setdefault(r, []).append(msgpack.unpackb(s, use_list=False))

    for r in d_n:
        d_n[r] = tuple(d_n[r])

    return d_n
