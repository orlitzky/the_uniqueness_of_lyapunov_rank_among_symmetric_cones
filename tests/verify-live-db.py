r"""
This script verifies that the live SQL database contains enough
information to verify the results in the paper. This is necessary
because the test suite is designed to pass regardless of how much
information is contained in the database (the tests should pass if the
database is completely empty).

This script can be executed directly,

    $ python verify-live-db.py

and it will return either success (0) or failure (1) after telling you
what it is doing.
"""

import msgpack
from random import randint
import sqlite3

from cones import SymmetricCone
from sql import *

if __name__ == "__main__":
    result = True

    mcd = max_cone_dim()
    mlrd = max_lorentz_rank_dim()

    def report():
        if result:
            print("ok", flush=True)
        else:
            print("fail", flush=True)


    print("Checking for the trivial cone in both tables... ", flush=True, end="")
    if not (have_cone_dim(0) and have_lorentz_rank_dim(0)):
        result = False
    report()

    print("All cone dimensions up to the max are present... ", flush=True, end="")
    for n in range(mcd + 1):
        if not have_cone_dim(n):
            result = False
    report()

    print("All lorentz rank dimensions up to the max are present... ", flush=True, end="")
    for n in range(mlrd + 1):
        if not have_lorentz_rank_dim(n):
            result = False
    report()

    print("Verifying the rank/dimension of some cones of a random size... ", flush=True, end="")
    if mcd < 0:
        result = False
    else:
        # decide on a dimension first; otherwise random ordering is
        # too slow.
        n = randint(0, mcd)
        conn = sqlite3.connect(LIVE_DATABASE)
        stmt = "SELECT rank,data FROM cones WHERE dim=? ORDER BY random() LIMIT 10000"
        with conn:
            rows = conn.execute(stmt, (n,)).fetchall()
        conn.close()

        for row in rows:
            K = SymmetricCone.deserialize(msgpack.unpackb(row[1], use_list=False))
            if not (K.dim == n and K.rank == row[0]):
                result = False
    report()

    print("Need max_cone_dim() >= 4 for Lemma 5... ", flush=True, end="")
    if mcd < 4:
        result = False
    report()

    print("Need max_cone_dim() >= 6 for non-Lorentz cones to exist... ", flush=True, end="")
    if mcd < 6:
        result = False
    report()

    print("Need max_cone_dim() >= 9 for Propositions 4 and 8... ", flush=True, end="")
    if mcd < 9:
        result = False
    report()

    print("Need max_cone_dim() >= 18 for Corollary 2... ", flush=True, end="")
    if mcd < 18:
        result = False
    report()

    print("Need max_cone_dim() >= 21 for lowerbound3b() and Theorem 5... ", flush=True, end="")
    if mcd < 21:
        result = False
    report()

    print("Need max_cone_dim() >= 27 for Example 1... ", flush=True, end="")
    if mcd < 27:
        result = False
    report()

    print("Need max_cone_dim() >= 36 for Proposition 9... ", flush=True, end="")
    if mcd < 36:
        result = False
    report()

    print("Need max_cone_dim() >= 39 for Proposition 7... ", flush=True, end="")
    if mcd < 39:
        result = False
    report()

    print("Need max_cone_dim() >= 40 for Corollary 3... ", flush=True, end="")
    if mcd < 40:
        result = False
    report()

    print("Need max_lorentz_rank_dim() >= 39 for Proposition 7... ", flush=True, end="")
    if mlrd < 39:
        result = False
    report()

    print("Need max_lorentz_rank_dim() >= 40 for Corollary 3... ", flush=True, end="")
    if mlrd < 40:
        result = False
    report()

    print("Confirming admissible_lorentz_ranks(60) directly... ", flush=True, end="")
    if mlrd < 60:
        result = False
    else:
        from partitions import _direct_lorentz_ranks
        actual = admissible_lorentz_ranks(60)
        expected = _direct_lorentz_ranks(60)
        if (set(actual) != set(expected)):
            result = False
    report()

    # This takes forever to get a random sample but there's no easy
    # way around it.
    print("There should be no hash collisions in the database... ", flush=True, end="")
    conn = sqlite3.connect(LIVE_DATABASE)
    expected = 1000000
    stmt = f"SELECT data FROM cones ORDER BY random() LIMIT {expected}"
    with conn:
        rows = conn.execute(stmt).fetchall()
    conn.close()

    # Dedupe via set(); if any were removed, there'll be fewer elements
    # than we started with.
    actual = len(set( hash(msgpack.unpackb(r[0], use_list=False))
                      for r in rows ))
    if not (actual == expected):
        result = False
    report()

    exit(int(not result))
