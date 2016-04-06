Installing and Basic Usage of XrayDB
=====================================

The X-ray database is held in the SQLite3 file ``xraydb.sqlite``.   For use
with SQLite, this file is all you need, as with::

   system~> sqlite3 xraydb.sqlite
   sqlite> .headers on
   sqlite> select * from elements where atomic_number=47;
   atomic_number|element|molar_mass|density
   47|Ag|107.868|10.48


Of course, the expectation is that you'd want to use this database within a
programming environment.   Currently, wrappers exist only for Python.
