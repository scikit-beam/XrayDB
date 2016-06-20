#!/usr/bin/env python
"""
SQLAlchemy wrapping of x-ray database for data from
     Elam et al, Chantler et al, Waasmaier and Kirfel

Main Class for full Database:  xrayDB
"""

import os
import time
import json
from collections import namedtuple
import numpy as np
from scipy.interpolate import interp1d, splrep, UnivariateSpline
from sqlalchemy import MetaData, create_engine
from sqlalchemy.orm import sessionmaker,  mapper, clear_mappers
from sqlalchemy.pool import SingletonThreadPool

# needed for py2exe?
import sqlalchemy.dialects.sqlite


XrayEdge = namedtuple('XrayEdge', ('edge', 'fyield', 'jump_ratio'))
XrayLine = namedtuple('XrayLine', ('energy', 'intensity', 'initial_level',
                           'final_level'))
CoreHoleWidth = namedtuple('CoreHoleWidth', ('Z', 'edge', 'width'))
ElementData = namedtuple('ElementData', ('Z', 'symbol', 'mass', 'density'))

def as_ndarray(obj):
    """make sure a float, int, list of floats or ints,
    or tuple of floats or ints, acts as a numpy array
    """
    if isinstance(obj, (float, int)):
        return np.array([obj])
    return np.asarray(obj)

def make_engine(dbname):
    return create_engine('sqlite:///%s' % (dbname),
                         poolclass=SingletonThreadPool)

def isxrayDB(dbname):
    """
    return whehter a file is a valid scan database:

    Parameters:
        dbname (string): name of Xray DB file

    Notes:
        must be a sqlite db file, with tables named 'elements',
        'photoabsorption', 'scattering', 'Coster_Kronig',
        'Chantler', 'Waasmaier', 'Version', and 'KeskiRahkonen_Krause'
    """
    _tables = ('Chantler', 'Waasmaier', 'Coster_Kronig',
               'KeskiRahkonen_Krause', 'Version',
               'elements', 'photoabsorption', 'scattering')
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
    "return json encoded value"
    if val is None or isinstance(val, (str, unicode)):
        return val
    return  json.dumps(val)


def elam_spline(xin, yin, yspl_in, x):
    """
    interpolate values from Elam photoabsorption and
    scattering tables, according to Elam, and following
    standard interpolation methods.  Calc borrowed from D. Dale.

    Parameters:
        xin (ndarray): x values for interpolation data
        yin (ndarray): y values for interpolation data
        yspl_in (ndarray): spline coefficients (second derivatives of y) for
                       interpolation data
        x (float or ndarray): x values to be evaluated at

    Returns:
        ndarray of interpolated values
    """
    x = as_ndarray(x)
    x[np.where(x < min(xin))] =  min(xin)
    x[np.where(x > max(xin))] =  max(xin)

    lo, hi = np.array([(np.flatnonzero(xin < e)[-1],
                        np.flatnonzero(xin > e)[0])
                       for e in x]).transpose()

    diff = xin[hi] - xin[lo]
    if any(diff <= 0):
        raise ValueError('x must be strictly increasing')
    a = (xin[hi] - x) / diff
    b = (x - xin[lo]) / diff
    return (a * yin[lo] + b * yin[hi] +
            (diff*diff/6) * ((a*a - 1) * a * yspl_in[lo] +
                             (b*b - 1) * b * yspl_in[hi] ))


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
    def __repr__(self):
        el = getattr(self, 'element', '??')
        line = getattr(self, 'siegbahn_symbol', '??')
        return "<%s(%s %s)>" % (self.__class__.__name__, el, line)

class WaasmaierTable(_BaseTable):
    (id, atomic_number, element, ion, offset, scale, exponents) = [None]*7
    def __repr__(self):
        el = getattr(self, 'ion', '??')
        return "<%s(%s)>" % (self.__class__.__name__, el)

class KeskiRahkonenKrauseTable(_BaseTable):
    (id, atomic_number, element, edge, width) = [None]*5
    def __repr__(self):
        el = getattr(self, 'element', '??')
        edge = getattr(self, 'edge', '??')
        return "<%s(%s %s)>" % (self.__class__.__name__, el, edge)

class ChantlerTable(_BaseTable):
    (id, element, sigma_mu, mue_f2, density,
     corr_henke, corr_cl35, corr_nucl,
     energy, f1, f2, mu_photo, mu_incoh, mu_total) = [None]*14

class XrayDB(object):
    "interface to Xray Data"
    def __init__(self, dbname='xraydb.sqlite', read_only=True):
        "connect to an existing database"
        if not os.path.exists(dbname):
            parent, child = os.path.split(__file__)
            dbname = os.path.join(parent, dbname)
            if not os.path.exists(dbname):
                raise IOError("Database '%s' not found!" % dbname)

        if not isxrayDB(dbname):
            raise ValueError("'%s' is not a valid X-ray Database file!" % dbname)

        self.dbname = dbname
        self.engine = make_engine(dbname)
        self.conn = self.engine.connect()
        kwargs = {}
        if read_only:
            kwargs = {'autoflush': True, 'autocommit':False}
            def readonly_flush(*args, **kwargs):
                return
            self.session = sessionmaker(bind=self.engine, **kwargs)()
            self.session.flush = readonly_flush
        else:
            self.session = sessionmaker(bind=self.engine, **kwargs)()

        self.metadata =  MetaData(self.engine)
        self.metadata.reflect()
        tables = self.tables = self.metadata.tables
        try:
            clear_mappers()
        except:
            pass
        mapper(ChantlerTable,            tables['Chantler'])
        mapper(WaasmaierTable,           tables['Waasmaier'])
        mapper(KeskiRahkonenKrauseTable, tables['KeskiRahkonen_Krause'])
        mapper(ElementsTable,            tables['elements'])
        mapper(XrayLevelsTable,          tables['xray_levels'])
        mapper(XrayTransitionsTable,     tables['xray_transitions'])
        mapper(CosterKronigTable,        tables['Coster_Kronig'])
        mapper(PhotoAbsorptionTable,     tables['photoabsorption'])
        mapper(ScatteringTable,          tables['scattering'])


    def close(self):
        "close session"
        self.session.flush()
        self.session.close()

    def query(self, *args, **kws):
        "generic query"
        return self.session.query(*args, **kws)

    def f0_ions(self, element=None):
        """
        return list of ion names supported for the .f0() calculation
        from Waasmaier and Kirfel

        Parameters:
            element (string, int, pr None):  atomic number, symbol, or ionic symbol
                    of scattering element.

        Returns:
            if element is None, all 211 ions are returned.
            if element is not None, the ions for that element are returned

        Example:
             >>> xdb = XrayDB()
             >>> xdb.f0_ions('Fe')
             ['Fe', 'Fe2+', 'Fe3+']

        Notes:
           Z values from 1 to 98 (and symbols 'H' to 'Cf') are supported.
        """
        rows = self.query(WaasmaierTable)
        if element is not None:
            rows = rows.filter(WaasmaierTable.element==self.symbol(element))
        return [str(r.ion) for r in rows.all()]

    def f0(self, ion, q):
        """
        return f0(q) -- elastic X-ray scattering factor from Waasmaier and Kirfel

        Parameters:
            ion (string, int, or None):  atomic number, symbol or ionic symbol
                  of scattering element.
            q (float, list, ndarray): value(s) of q for scattering factors

        Returns:
            ndarray of elastic scattering factors

        Notes:
            q = sin(theta) / lambda, where theta = incident angle, and
            lambda = X-ray wavelength
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
            q = as_ndarray(q)
            f0 = row.offset
            for s, e in zip(json.loads(row.scale), json.loads(row.exponents)):
                f0 += s * np.exp(-e*q*q)
            return f0

    def _getChantler(self, element, energy, column='f1', smoothing=1):
        """
        return energy-dependent data from Chantler table

        Parameters:
            element (string or int): atomic number or symbol.
            eneregy (float or ndarray):
        columns: f1, f2, mu_photo, mu_incoh, mu_total

        Notes:
           this function is meant for internal use.
        """
        tab = ChantlerTable
        row = self.query(tab)
        row = row.filter(tab.element==self.symbol(element)).all()
        if len(row) > 0:
            row = row[0]
        if isinstance(row, tab):
            energy = as_ndarray(energy)
            emin, emax = min(energy), max(energy)
            # te = self.chantler_energies(element, emin=emin, emax=emax)
            te = np.array(json.loads(row.energy))
            nemin = max(0, -5 + max(np.where(te<=emin)[0]))
            nemax = min(len(te), 6 + max(np.where(te<=emax)[0]))
            region = np.arange(nemin, nemax)
            te = te[region]
            if column == 'mu':
                column = 'mu_total'
            ty = np.array(json.loads(getattr(row, column)))[region]
            if column == 'f1':
                out = UnivariateSpline(te, ty, s=smoothing)(energy)
            else:
                out = np.exp(np.interp(np.log(energy),
                                       np.log(te),
                                       np.log(ty)))
            if isinstance(out, np.ndarray) and len(out) == 1:
                return out[0]
            return out

    def chantler_energies(self, element, emin=0, emax=1.e9):
        """
        return array of energies (in eV) at which data is
        tabulated in the Chantler tables for a particular element.

        Parameters:
            element (string or int): atomic number or symbol
            emin (float): minimum energy (in eV) [0]
            emax (float): maximum energy (in eV) [1.e9]

        Returns:
            ndarray of energies
        """
        tab = ChantlerTable
        row = self.query(tab).filter(tab.element==self.symbol(element)).all()
        if len(row) > 0:
            row = row[0]
        if not isinstance(row, tab):
            return None
        te = np.array(json.loads(row.energy))
        tf1 = np.array(json.loads(row.f1))
        tf2 = np.array(json.loads(row.f2))

        if emin <= min(te):
            nemin = 0
        else:
            nemin = max(0,  -2 + max(np.where(te<=emin)[0]))
        if emax > max(te):
            nemax = len(te)
        else:
            nemax = min(len(te), 2 + max(np.where(te<=emax)[0]))
        region = np.arange(nemin, nemax)
        return te[region] # , tf1[region], tf2[region]

    def f1_chantler(self, element, energy, **kws):
        """
        returns f1 -- real part of anomalous X-ray scattering factor
        for selected input energy (or energies) in eV.

        Parameters:
            element (string or int): atomic number or symbol
            energy (float or ndarray): energies (in eV).

        Returns:
            ndarray of anomalous scattering factor
        """
        return self._getChantler(element, energy, column='f1', **kws)

    def f2_chantler(self, element, energy, **kws):
        """
        returns f2 -- imaginary part of anomalous X-ray scattering factor
        for selected input energy (or energies) in eV.

        Parameters:
            element (string or int): atomic number or symbol
            energy (float or ndarray): energies (in eV).

        Returns:
            ndarray of anomalous scattering factor
        """
        return self._getChantler(element, energy, column='f2', **kws)

    def mu_chantler(self, element, energy, incoh=False, photo=False):
        """
        returns X-ray mass attenuation coefficient, mu/rho in cm^2/gr
        for selected input energy (or energies) in eV.
        default is to return total attenuation coefficient.

        Parameters:
            element (string or int): atomic number or symbol
            energy (float or ndarray): energies (in eV).
            photo (bool): return only the photo-electric contribution [False]
            incoh (bool): return only the incoherent contribution [False]

        Returns:
            ndarray of mass attenuation coefficient.
        """
        col = 'mu_total'
        if photo:
            col = 'mu_photo'
        elif incoh:
            col = 'mu_incoh'
        return self._getChantler(element, energy, column=col)

    def _getElementData(self, element):
        "return data from elements table: internal use"
        elem = ElementsTable.element
        if isinstance(element, int):
            elem = ElementsTable.atomic_number
        row = self.query(ElementsTable).filter(elem==element).all()
        if len(row) > 0:
            row = row[0]
        return ElementData(int(row.atomic_number),
                           row.element.title(),
                           row.molar_mass, row.density)

    def zofsym(self, element):
        """
        return element's atomic number

        Parameters:
            element (string or int): atomic number or symbol

        Returns:
            atomic number (integer)
        """
        return self._getElementData(element).Z

    def symbol(self, element):
        """
        return element symbol

        Parameters:
            element (string or int): atomic number or symbol

        Returns:
            element symbol (string)
        """
        return self._getElementData(element).symbol

    def molar_mass(self, element):
        """
        return molar mass of element

        Parameters:
            element (string or int): atomic number or symbol

        Returns:
            molar mass of element (float) in amu
        """
        return self._getElementData(element).mass

    def density(self, element):
        """
        return density of pure element

        Parameters:
            element (string or int): atomic number or symbol

        Returns:
            density of element (float) in grams/cm^3
       """
        return self._getElementData(element).density

    def xray_edges(self, element):
        """
        returns dictionary of X-ray absorption edge energy (in eV),
        fluorescence yield, and jump ratio for an element.

        Parameters:
            element (string or int): atomic number or symbol

        Returns:
            dictionary with keys of edge (iupac symbol), and values of
            XrayEdge (namedtuple of (energy, fyield, edge_jump))
        """
        element = self.symbol(element)
        tab = XrayLevelsTable
        out = {}
        for r in self.query(tab).filter(tab.element==element).all():
            out[str(r.iupac_symbol)] = XrayEdge(r.absorption_edge,
                                                r.fluorescence_yield,
                                                r.jump_ratio)
        return out

    def xray_edge(self, element, edge):
        """
        returns XrayEdge for an element and edge

        Parameters:
            element (string or int): atomic number or symbol
            edge (string):  X-ray edge

        Returns:
            XrayEdge (namedtuple of (energy, fyield, edge_jump))
        """
        edges = self.xray_edges(element)
        edge = edge.title()
        if edge in edges:
            return edges[edge]

    def xray_lines(self, element, initial_level=None, excitation_energy=None):
        """
        returns dictionary of X-ray emission lines of an element, with

        Parameters:
            initial_level (string or list/tuple of string):  initial level(s) to
                 limit output.
            excitation_energy (float): energy of excitation, limit output those
<                 excited by X-rays of this energy (in eV).

        Returns:
            dict with keys of lines (iupac symbol) and values of Xray Lines

        Note:
            if both excitation_energy and initial_level are given, excitation_level
            will limit output
        """
        element = self.symbol(element)
        tab = XrayTransitionsTable
        row = self.query(tab).filter(tab.element==element)
        if excitation_energy is not None:
            initial_level = []
            for ilevel, dat in self.xray_edges(element).items():
                if dat[0] < excitation_energy:
                    initial_level.append(ilevel.title())

        if initial_level is not None:
            if isinstance(initial_level, (list, tuple)):
                row = row.filter(tab.initial_level.in_(initial_level))
            else:
                row = row.filter(tab.initial_level==initial_level.title())
        out = {}
        for r in row.all():
            out[str(r.siegbahn_symbol)] = XrayLine(r.emission_energy, r.intensity,
                                               r.initial_level, r.final_level)
        return out

    def xray_line_strengths(self, element, excitation_energy=None):
        """
        return the absolute line strength in cm^2/gr for all available lines

        Parameters:
            element (string or int): Atomic symbol or number for element
            excitation_energy (float): incident energy, in eV

        Returns:
            dict of elemental line with fluorescence cross section in cm2/gr.

        """
        out = {}
        for label, eline in self.xray_lines(element, excitation_energy=excitation_energy).items():
            edge = self.xray_edge(element, eline.initial_level)
            if edge is None and ',' in eline.initial_level:
                ilevel, extra = eline.initial_level.split(',')
                edge = self.xray_edge(element, ilevel)
            if edge is not None:
                mu = self.mu_elam(element, [edge.absorption_edge*(0.999),
                                            edge.absorption_edge*(1.001)],
                                  kind='photo')
                out[label] = (mu[1]-mu[0]) * eline.intensity * edge.fluorescence_yield
        return out

    def CK_probability(self, element, initial, final, total=True):
        """
        return Coster-Kronig transition probability for an element and
        initial/final levels

        Parameters:
            element (string or int): Atomic symbol or number for element
            initial (string):  initial level
            final (string):  final level
            total (bool): whether to return total or partial probability

        Returns:
            float transition probability

        Example:
            >>> xdb = XrayDB()
            >>> xdb.CK_probability('Cu', 'L1', 'L3', total=True)
            0.681
        """
        element = self.symbol(element)
        tab = CosterKronigTable
        row = self.query(tab).filter(
            tab.element==element).filter(
            tab.initial_level==initial.title()).filter(
            tab.final_level==final.title()).all()
        if len(row) > 0:
            row = row[0]
        if isinstance(row, tab):
            if total:
                return row.total_transition_probability
            else:
                return row.transition_probability

    def corehole_width(self, element=None, edge=None):
        """
        returns core hole width for an element and edge

        Parameters:
         element (string, integer, or None): atomic number or symbol for element
         edge (sring or None): edge for hole.

        Returns:
          list of named tuples (Z, edge, width) for element and edge
          with widths in eV.

        Notes:
          if edge is None, values for all edges for the element are returned.
          if element is None, values for all elements are returned.
        """
        tab = KeskiRahkonenKrauseTable
        rows = self.query(tab)
        if element is not None:
            rows = rows.filter(tab.element==self.symbol(element))
        if edge is not None:
            rows = rows.filter(tab.edge==edge.title())
        return [CoreHoleWidth(r.atomic_number, r.edge, r.width)
                for r in rows.all()]

    def Elam_CrossSection(self, element, energies, kind='photo'):
        """
        returns Elam Cross Section values for an element and energies

        Parameters:
            element (string or int):  atomic number or symbol for element
            energies (float or ndarray): energies (in eV) to calculate cross-sections
            kind (string):  one of 'photo', 'coh', and 'incoh' for photo-absorption,
                  coherent scattering, and incoherent scattering cross sections,
                  respectively. Default is 'photo'.

        Returns:
            ndarray of scattering data

        Data from Elam, Ravel, and Sieber.
        """
        element = self.symbol(element)
        energies = 1.0 * as_ndarray(energies)

        kind = kind.lower()
        if kind not in ('coh', 'incoh', 'photo'):
            raise ValueError('unknown cross section kind=%s' % kind)

        tab = ScatteringTable
        if kind == 'photo':
            tab = PhotoAbsorptionTable

        row = self.query(tab).filter(tab.element==element).all()
        if len(row) > 0:
            row = row[0]
        if not isinstance(row, tab):
            return None
        tab_lne = np.array(json.loads(row.log_energy))
        if kind.startswith('coh'):
            tab_val = np.array(json.loads(row.log_coherent_scatter))
            tab_spl = np.array(json.loads(row.log_coherent_scatter_spline))
        elif kind.startswith('incoh'):
            tab_val = np.array(json.loads(row.log_incoherent_scatter))
            tab_spl = np.array(json.loads(row.log_incoherent_scatter_spline))
        else:
            tab_val = np.array(json.loads(row.log_photoabsorption))
            tab_spl = np.array(json.loads(row.log_photoabsorption_spline))

        emin_tab = 10*int(0.102*np.exp(tab_lne[0]))
        energies[np.where(energies < emin_tab)] = emin_tab
        out = np.exp(elam_spline(tab_lne, tab_val, tab_spl, np.log(energies)))
        if len(out) == 1:
            return out[0]
        return out

    def mu_elam(self, element, energies, kind='total'):
        """
        returns attenuation cross section for an element at energies (in eV)

        Parameters:
            element (string or int):  atomic number or symbol for element
            energies (float or ndarray): energies (in eV) to calculate cross-sections
            kind (string):  one of 'photo' or 'total' for photo-electric or
                  total attenuation, respectively.  Default is 'total'.

        Returns:
           ndarray of scattering values in units of cm^2 / gr

        Data from Elam, Ravel, and Sieber.
        """
        calc = self.Elam_CrossSection
        xsec = calc(element, energies, kind='photo')
        if kind.lower().startswith('tot'):
            xsec += calc(element, energies, kind='coh')
            xsec += calc(element, energies, kind='incoh')
        return xsec

    def coherent_cross_section_elam(self, element, energies):
        """
        returns coherenet scattering cross section for an element
        at energies (in eV)

        Parameters:
            element (string or int):  atomic number or symbol for element
            energies (float or ndarray): energies (in eV) to calculate cross-sections

        Returns:
           values in units of cm^2 / gr

        Data from Elam, Ravel, and Sieber.
        """
        return self.Elam_CrossSection(element, energies, kind='coh')

    def incoherent_cross_section_elam(self, element, energies):
        """
        returns incoherenet scattering cross section for an element
        at energies (in eV)


        Parameters:
            element (string or int):  atomic number or symbol for element
            energies (float or ndarray): energies (in eV) to calculate cross-sections

        Returns:
           values in units of cm^2 / gr

        Data from Elam, Ravel, and Sieber.
        """
        return self.Elam_CrossSection(element, energies, kind='incoh')
