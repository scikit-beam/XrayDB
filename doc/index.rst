.. xraydb documentation master file

X-ray DB: X-ray Reference Data in SQLite
========================================

XrayDB provides Atomic data, characteristic X-ray energies, and X-ray cross
sections for the elements in an SQLite3 database, ``xraydb.sqlite``.  This
file can be used directly with SQLite using standard SQL, or or with any of
several programming language that has an SQLite library.  See
http://sqlite.org for further details on SQLite.  As of this writing, a
wrapper for Python (http://python.org) is provided.

Because some of the components of the database hold arrays of numbers
(for example, coefficients for interpolation), the arrays are stored in the
database as JSON-encoded strings, and will need to be unpacked to be used.
For reference, see http://json.org.

The project began with the data from the compilation of basic atomic
properties and X-ray absorption edge energies, emission energies, and
absorption cross sections from :cite:ts:`Elam_etal`.  More data has since
been added from several sources.  Energy widths of core holes for excited
electronic levels from :cite:ts:`Keski_Krause` is included.  Elastic X-ray
scattering data, :math:`f_0(q)` is taken from :cite:ts:`Waasmaier_Kirfel`.
Resonant scattering cross sections :math:`f'(E)` and :math:`f''(E)` and and
absorption cross sections from :cite:ts:`Chantler` (as from
http://www.nist.gov/pml/data/ffast/index.cfm) are also included.

.. toctree::
   :maxdepth: 2

   installation
   overview
   dbschema
   python
