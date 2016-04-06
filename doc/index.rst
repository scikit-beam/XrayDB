.. xraydb documentation master file

X-ray DB: X-ray Reference Data in SQLite
========================================

X-ray DB provides ``xraydb.sqlite``, an SQLite3 database of atomic data,
characteristic X-ray energies and X-ray cross sections for the elements.
The project began with the data from the compilation of Elam, Ravel, and
Sieber.  More data has since been added.

The xraydb.sqlite can be used directly with sqlite or with any programming
language that has an SQLite library.  See http://sqlite.org for further details.

Because some of the components of the database hold arrays of numbers
(for example, coefficients for interpolation), the arrays are stored in the
database as JSON-encoded strings, and will need to be unpacked to be used.
For reference, see http://json.org.

For some programming languages, higher-level wrappers may be provided to give an
even easier interface to these data.

.. toctree::
   :maxdepth: 2

   installation
   python
