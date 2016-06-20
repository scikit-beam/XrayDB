Using XrayDB from Python
=========================


The wrappers/python directly contains a Python module for xraydb.  This
module gives a higher-level wrapping of the XrayDB, including the
conversion of data from json-encoded data to numpy arrays.  The module
requires the json, numpy, and sqlalchemy modules, and can be installed
with::

    python setup.py install



Using the XrayDB database
============================

.. module:: xraydb

To use the XrayDB from python, create an instance, and start using it:

    >>> from xraydb import XrayDB
    >>> xdb = XrayDB()
    >>> xdb.xray_edge('Ag', 'K')
    XrayEdge(edge=25514.0, fyield=0.821892, jump_ratio=6.334)



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
