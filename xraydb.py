#!/usr/bin/env python
"""
SQLAlchemy wrapping of x-ray database for Elam and Chantler data

Main Class for full Database:  XrayDB
"""

import os
import json
import numpy as np
import time

from sqlalchemy import MetaData, and_, create_engine, \
     Table, Column, Integer, Float, String, Text, DateTime, ForeignKey

from sqlalchemy.orm import sessionmaker,  mapper, clear_mappers, relationship
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm.exc import  NoResultFound
from sqlalchemy.pool import SingletonThreadPool

# needed for py2exe?
import sqlalchemy.dialects.sqlite

def make_engine(dbname):
    return create_engine('sqlite:///%s' % (dbname),
                         poolclass=SingletonThreadPool)

def isxrayDB(dbname):
    """test if a file is a valid scan database:
    must be a sqlite db file, with tables named
    'Coster_Kronig', 'elements', 'photoabsorption', 'scattering'

    """
    _tables = ('Chantler', 'Waasmaier', 'Coster_Kronig', 'elements',
               'photoabsorption', 'scattering')
    result = False
    try:
        engine = make_engine(dbname)
        meta = MetaData(engine)
        meta.reflect()
        result = all([t in meta.tables for t in _tables])
    except:
        pass
    return result

def json_encode(val):
    "simple wrapper around json.dumps"
    if val is None or isinstance(val, (str, unicode)):
        return val
    return  json.dumps(val)

class DBException(Exception):
    """DB Access Exception: General Errors"""
    def __init__(self, msg):
        Exception.__init__(self)
        self.msg = msg
    def __str__(self):
        return self.msg

class _BaseTable(object):
    "generic class to encapsulate SQLAlchemy table"
    def __repr__(self):
        el = getattr(self, 'element', '??')
        return "<%s(%s)>" % (self.__class__.__name__, el)

class CosterKronigTable(_BaseTable):
    "positioners table"
    (id, element, initial_level, final_level,
     transition_probability, total_transition_probability) = [None]*6

class ElementsTable(_BaseTable):
    (atomic_number, element, molar_mass, density) = [None]*4

class PhotoAbsorptionTable(_BaseTable):
    (id, element, log_energy,
     log_photoabsorption, log_photoabsorption_spline) = [None]*5

class ScatteringTable(_BaseTable):
    (id, element, log_energy,
     log_coherent_scatter, log_coherent_scatter_spline,
     log_incoherent_scatter, log_incoherent_scatter_spline) = [None]*7

class XrayLevelsTable(_BaseTable):
    (id, element,  iupac_symbol,
     absorption_edge, fluorescence_yield, jump_ratio) = [None]*6
    def __repr__(self):
        el = getattr(self, 'element', '??')
        edge= getattr(self, 'iupac_symbol', '??')
        return "<%s(%s %s)>" % (self.__class__.__name__, el, edge)


class XrayTransitionsTable(_BaseTable):
    (id, element, iupac_symbol, siegbahn_symbol, initial_level,
     final_level, emission_energy, intensity) = [None]*8

class WaasmaierTable(_BaseTable):
    (id, atomic_number, element, ion, offset, scale, exponents) = [None]*7
    def __repr__(self):
        el = getattr(self, 'ion', '??')
        return "<%s(%s)>" % (self.__class__.__name__, el)

class ChantlerTable(_BaseTable):
    (id, element, sigma_mu, mue_f2, density,
     energy, f1, f2, mu_photo, mu_incoh, mu_total) = [None]*11

class xrayDB(object):
    "interface to Xray Data"
    def __init__(self, dbname='xrayref.db'):
        "connect to an existing database"
        if not os.path.exists(dbname):
            raise IOError("Database '%s' not found!" % dbname)

        if not isxrayDB(dbname):
            raise ValueError("'%s' is not a valid X-ray Database file!" % dbname)

        self.dbname = dbname
        self.engine = make_engine(dbname)
        self.conn = self.engine.connect()
        self.session = sessionmaker(bind=self.engine)()
        self.metadata =  MetaData(self.engine)
        self.metadata.reflect()
        tables = self.tables = self.metadata.tables
        try:
            clear_mappers()
        except:
            pass
        mapper(ChantlerTable,        tables['Chantler'])
        mapper(WaasmaierTable,       tables['Waasmaier'])
        mapper(CosterKronigTable,    tables['Coster_Kronig'])
        mapper(ElementsTable,        tables['elements'])
        mapper(PhotoAbsorptionTable, tables['photoabsorption'])
        mapper(ScatteringTable,      tables['scattering'])
        mapper(XrayLevelsTable,      tables['xray_levels'])
        mapper(XrayTransitionsTable, tables['xray_transitions'])

    def close(self):
        "close session"
        self.session.flush()
        self.session.close()

    def query(self, *args, **kws):
        "generic query"
        return self.session.query(*args, **kws)

    def f0_ions(self, element=None):
        """return list of ion names supported for the .f0() calculation
        from Waasmaier and Kirfel

        if element is None, all 211 ions are returned.

        If element is not None, the ions for that element (atomic symbol) are returned
        """
        rows = self.query(WaasmaierTable)
        if element is not None:
            if isinstance(element, int):
                rows = rows.filter(WaasmaierTable.atomic_number==element)
            else:
                rows = rows.filter(WaasmaierTable.element==element.title())
        return [str(r.ion) for r in rows.all()]

    def f0(self, q, ion):
        """Calculate f0(q) -- elastic x-ray scattering factor
        from Waasmaier and Kirfel

        arguments
        ---------
        q: single q value, list, tuple, or numpy array of q value
             q = sin(theta) / lambda
             theta = incident angle, lambda = x-ray wavelength

        ion:  atomic number, atomic symbol or ionic symbol
              (case insensitive) of scatterer

        Z values from 1 to 98 (and symbols 'H' to 'Cf') are supported.
        The list of ionic symbols can be read with the function .f0_ions()
        """
        tab = WaasmaierTable
        row = self.query(tab)
        if isinstance(ion, int):
            row = row.filter(tab.atomic_number==ion).all()
        else:
            row = row.filter(tab.ion==ion.title()).all()
        if len(row) > 0:
            row = row[0]
        if isinstance(row, tab):
            if isinstance(q, (tuple, list)):
                q = np.array(q)
            f0 = row.offset
            for s, e in zip(json.loads(row.scale), json.loads(row.exponents)):
                f0 += s * np.exp(-e*q*q)
            return f0

    def _getChantler(self, energy, element, column='f1'):
        """return energy-dependent data from Chantler table
        columns: f1, f2, mu_photo, mu_incoh, mu_total
        """
        tab = ChantlerTable
        row = self.query(tab)
        if isinstance(element, int):
            row = row.filter(tab.id==element).all()
        else:
            row = row.filter(tab.element==element.title()).all()
        if len(row) > 0:
            row = row[0]
        if isinstance(row, tab):
            if isinstance(energy, (tuple, list)):
                energy = np.array(energy)
            if column == 'mu':
                column = 'mu_total'
            if hasattr(row, column):
                return np.interp(energy, np.array(json.loads(row.energy)),
                                 np.array(json.loads(getattr(row, column))))

    def f1(self, energy, element):
        """returns f1 -- real part of anomalous x-ray scattering factor
        for selected input energy (or energies) in eV.
        """
        return self._getChantler(energy, element, column='f1')

    def f2(self, energy, element):
        """returns f2 -- imaginary part of anomalous x-ray scattering factor
        for selected input energy (or energies) in eV.
        """
        return self._getChantler(energy, element, column='f2')

    def mu(self, energy, element, incoh=False, photo=False):
        """returns mu/rho in cm^2/gr -- x-ray mass attenuation coefficient
        for selected input energy (or energies) in eV.
        default is to return total attenuation coefficient.
        use
          photo=True to return only the photo-electric contribution or
          incoh=True to return on the incoherent contribution
        """
        col = 'mu_total'
        if incoh:
            col = 'mu_incoh'
        elif photo:
            col = 'mu_photo'
        return self._getChantler(energy, element, column=col)

    def _getElementData(self, element):
        "get data from elements table"
        tab = ElementsTable
        row = self.query(tab)
        if isinstance(element, int):
            row = row.filter(tab.atomic_number==element).all()
        else:
            row = row.filter(tab.element==element.title()).all()
        if len(row) > 0: row = row[0]
        return row

    def zofsym(self, element):
        "return z for element name"
        return int(self._getElementData(element).atomic_number)

    def symbol(self, z):
        "return element symbol from z"
        return self._getElementData(z).element

    def molar_mass(self, element):
        "return molar mass of element"
        return self._getElementData(element).molar_mass

    def density(self, element):
        "return density of pure element"
        return self._getElementData(element).density


    def xray_edges(self, element):
        """returns dictionary of all x-ray absorption edge energy (in eV), fluorescence yield, and jump ratio
        for an element dictionary has keys of edge (iupac symol), each containing a tuple of
               (energy, fluorescence_yield, edge_jump)
        """
        if isinstance(element, int):
            element = self.symbol(element)
        tab = XrayLevelsTable
        out = {}
        for r in self.query(tab).filter(tab.element==element.title()).all():
            out[str(r.iupac_symbol)] = (r.absorption_edge,
                                        r.fluorescence_yield,
                                        r.jump_ratio)
        return out

    def xray_edge(self, element, edge):
        """returns tuple of
               (energy, fluorescence_yield, edge_jump)
        for an x-ray absorption edge
        """
        edges = self.xray_edges(element)
        if edge in edges:
            return edges[edge]

    def xray_lines(self, element): # iupac=None, siegbahn=None,   excitation_energy=None):
        """returns x-ray emission line energy (in eV), line intensity, iupac symbol,
           siegbahn symbol, initial level and final level
        for an element
        """
        if isinstance(element, int):
            element = self.symbol(element)
        tab = XrayTransitionsTable
#         if iupac is not None:
#             row = row.filter(tab.iupac_symbol==iupac)
#         if siegbahn is not None:
#             row = row.filter(tab.siegbahn_symbol==siegbahn)
        out = {}
        for r in self.query(tab).filter(tab.element==element.title()).all():
            out[str(r.siegbahn_symbol)] = (r.emission_energy,
                                        r.intensity,
                                        r.initial_level,
                                        r.final_level, r.iupac_symbol)
        return out
