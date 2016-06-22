Using XrayDB from Python
=========================


The wrappers/python directly contains a Python module for xraydb.  This
module gives a higher-level wrapping of the XrayDB, including the
conversion of data from json-encoded data to numpy arrays.  The module
requires the json, numpy, and sqlalchemy modules, and can be installed
with::

    python setup.py install


XrayDB module
---------------

.. module:: xraydb

To use the XrayDB from python, create an instance, and start using it:

    >>> from xraydb import XrayDB
    >>> xdb = XrayDB()
    >>> xdb.xray_edge('Ag', 'K')
    XrayEdge(edge=25514.0, fyield=0.821892, jump_ratio=6.334)


.. index:: XrayDB methods
.. _xraydb-methods_table:

    Table of XrayDB methods for Atomic and X-ray data for the elements.
    calculate and return some element-specific properties, given the
    element symbol or atomic number.  Most data extends to Z=98 (Cf), but
    much data for elements with atomic number > 92 (U) may not be
    available, and may not be very reliable when provided.  Except where
    noted, the data comes :cite:ts:`Elam_etal`.

     ========================== =============================================================
      function                    description
     ========================== =============================================================
      :meth:`atomic_number`      atomic number from symbol
      :meth:`atomic_symbol`      atomic symbol from number
      :meth:`atomic_mass`        atomic mass
      :meth:`atomic_density`     atomic density (for pure element)
      :meth:`xray_edge`          xray edge data for a particular element and edge
      :meth:`xray_line`          xray emission line data for an element and line
      :meth:`xray_edges`         dictionary of all X-ray edges data for an element
      :meth:`xray_lines`         dictionary of all X-ray emission line data for an element
      :meth:`fluo_yield`         fluorescence yield and weighted line energy
      :meth:`core_width`         core level width for an element and edge (:cite:ts:`Keski_Krause`)
      :meth:`mu_elam`            absorption cross section
      :meth:`coherent_xsec`      coherent cross section
      :meth:`incoherent_xsec`    incoherent cross section
      :meth:`f0`                 elastic scattering factor (:cite:ts:`Waasmaier_Kirfel`)
      :meth:`f0_ions`            list of valid "ions" for :meth:`f0`  (:cite:ts:`Waasmaier_Kirfel`)
      :meth:`chantler_energies`  energies of tabulation for Chantler data (:cite:ts:`Chantler`)
      :meth:`f1_chantler`        f'  anomalous factor (:cite:ts:`Chantler`)
      :meth:`f2_chantler`        f'' anomalous factor (:cite:ts:`Chantler`)
      :meth:`mu_chantler`        absorption cross section (:cite:ts:`Chantler`)
     ========================== =============================================================


.. autoclass:: XrayDB

    .. automethod:: zofsym

    .. automethod:: symbol

    .. automethod:: molar_mass

    .. automethod:: density

    .. automethod:: xray_edges

    .. automethod:: xray_edge

    .. automethod:: xray_lines

    .. automethod:: xray_line_strengths

    .. automethod:: CK_probability

    .. automethod:: corehole_width

    .. automethod:: Elam_CrossSection

    .. automethod:: mu_elam

    .. automethod:: coherent_cross_section_elam

    .. automethod:: incoherent_cross_section_elam

    .. automethod:: f0_ions

    .. automethod:: f0

    .. automethod:: chantler_energies

    .. automethod:: f1_chantler

    .. automethod:: f2_chantler

    .. automethod:: mu_chantler
