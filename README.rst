
XrayDB: X-ray Reference Data in SQLite
=======================================

See http://scikit-beam.github.io/XrayDB/

XrayDB provides atomic data, characteristic X-ray energies, and X-ray cross
sections for the elements in an SQLite3 database file ``xraydb.sqlite``.
This file can be used directly with sqlite3 using standard SQL or with any
programming language that has an interface to SQLite.  See
http://sqlite.org for further details on SQLite.  A Python module providing
an interface is also provided.

Because some of the components of the database hold arrays of numbers (for
example, coefficients for interpolation), the arrays are stored in the
database as JSON-encoded strings, and will need to be unpacked to be used.
For reference, see http://json.org.

The project began with the data from the compilation of basic atomic
properties and X-ray absorption edge energies, emission energies, and
absorption cross sections from Elam, Ravel, and Sieber[1], who assembled
data from a several sources.  More data has since been added additional
sources.  Energy widths of core holes for excited electronic levels from
Keski-Rahkonen and Krause[2] is included.  Elastic X-ray scattering data,
f0(q) is taken from Waasmaier and Kirfel[3].  Resonant scattering cross
sections f'(E) and f''(E) and absorption cross sections from Chantler[4]
(as from http://www.nist.gov/pml/data/ffast/index.cfm) are also included.

In general, cross sections are in cm^2/gr, and energies are given in eV.
Energy-dependent data for cross sections is typically valid in the range
from about 250 eV to about 200,000 eV.

See the documentation for the database schema, and for further references.

References
-----------

The data from Elam, Ravel and Sieber is itself a compilation from other
sources. The Data/elam.dat file (and the documentation above) for a more
complete set of references.

[1] W. T. Elam, B. D. Ravel and J. R. Sieber, Radiation Physics and Chemistry 63 (2),
    pp121–128 (2002) [http://dx.doi.org/10.1016/S0969-806X(01)00227-4].
[2] O. Keski-Rahkonen and M. O. Krause, Atomic Data and Nuclear Data Tables
    14, pp139-146 (1974) [http://dx.doi.org/10.1016/S0092-640X(74)80020-3]
[3] D. Waasmaier and A. Kirfel, Acta Crystallographica A51, pp416-431 (1995)
    [http://dx.doi.org/10.1107/S0108767394013292].
[4] C. T. Chantler, Journal of Physical and Chemical Reference Data 29,
    pp597-1048 (2000) [http://dx.doi.org/10.1063/1.1321055]
