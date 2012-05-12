#!/usr/bin/env python
"""
SQLAlchemy wrapping of x-ray database for Elam and Chantler data

Main Class for full Database:  XrayDB
"""

import os
import json
import epics
import time
import socket

from sqlalchemy import MetaData, and_, create_engine, \
     Table, Column, Integer, Float, String, Text, DateTime, ForeignKey

from sqlalchemy.orm import sessionmaker,  mapper, clear_mappers, relationship
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm.exc import  NoResultFound
from sqlalchemy.pool import SingletonThreadPool

# needed for py2exe?
import sqlalchemy.dialects.sqlite

def make_engine(dbname, server):
    return create_engine('sqlite:///%s' % (dbname),
                         poolclass=SingletonThreadPool)

def isXrayDB(dbname, server='sqlite'):
    """test if a file is a valid scan database:
    must be a sqlite db file, with tables named
       'postioners', 'detectors', and 'scans'
    """
    _tables = ('Coster_Kronig', 'elements', 'photoabsorption', 'scattering')
    result = False
    try:
        engine = make_engine(dbname, server)
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

def None_or_one(val, msg='Expected 1 or None result'):
    """expect result (as from query.all() to return
    either None or exactly one result
    """
    if len(val) == 1:
        return val[0]
    elif len(val) == 0:
        return None
    else:
        raise DBException(msg)


class DBException(Exception):
    """DB Access Exception: General Errors"""
    def __init__(self, msg):
        Exception.__init__(self)
        self.msg = msg
    def __str__(self):
        return self.msg

def StrCol(name, size=None, **kws):
    val = Text
    if size is not None: val = String(size)
    return Column(name, val, **kws)

def FloatCol(name, **kws):
    return Column(name, Float, **kws)

def NamedTable(tablename, metadata, keyid='id',
               nameid='element', name=True, cols=None):
    args  = [Column(keyid, Integer, primary_key=True)]
    if name:
        args.append(StrCol(nameid, size=16, nullable=False, unique=True))
    if cols is not None:
        args.extend(cols)
    return Table(tablename, metadata, *args)


def create_newdb(dbname, server='sqlite'):
    engine  = make_engine(dbname, server)
    metadata =  MetaData(engine)
    # print dbname, engine, metadataa
    ck = NamedTable('Coster_Kroning', metadata,
                    cols=[StrCol('initial_level'),
                          StrCol('final_level'),
                          FloatCol('transition_probability'),
                          FloatCold('total_transition_probability')])

    el = NamedTable('elements', metadata, keyid='atomic_number',
                    cols=[FloatCol('molar_mass'), FloatCol('density')])

    pa = NamedTable('photoabsorption', metadata,
                    cols=[StrCol('log_energy'), StrCol('log_photoabsorption'),
                          StrCol('log_photoabsorption_spline')])

    sc = NamedTable('scattering', metadata,
                    cols=[StrCol('log_energy'),
                          StrCol('log_coherent_scatter'),
                          StrCol('log_coherent_scatter_spline'),
                          StrCol('log_incoherent_scatter'),
                          StrCol('log_incoherent_scatter_spline')])

    xl = NamedTable('xray_levels', metadata,
                    cols=[StrCol('iupac_symbol'),
                          FloatCol('absorption_edge'),
                          FloatCol('fluorescence_yield'),
                          FloatCol('jump_ratio')])
    xt = NamedTable('xray_transitions', metadata,
                    cols=[StrCol('iupac_symbol'),
                          StrCol('siegbahn_symbol'),
                          StrCol('initial_level'),
                          StrCol('final_level'),
                          FloatCol('emission_energy'),
                          FloatCol('intensity')])


    metadata.create_all()
    session = sessionmaker(bind=engine)()
    session.commit()

class _BaseTable(object):
    "generic class to encapsulate SQLAlchemy table"
    def __repr__(self):
        el = getattr(self, 'element', '??')
        return "<%s(%s)>" % (self.__class__.__name__, el)

class CKTable(_BaseTable):
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
    
class XrayTransitionsTable(_BaseTable):
    (id, element, iupac_symbol, siegbahn_symbol, initial_level,
     final_level, emission_energy, intensity) = [None]*8

class WaasmaierTable(_BaseTable):
    (id, atomic_number, element, ion, offset, scale, exponents) = [None]*7
    
class ChantlerTable(_BaseTable):
    (id, element, sigma_mu, mue_f2, density,
     energy, f1, f2, mu_photo, mu_incoh, mu_total) = [None]*11

class xrayDB(object):
    "interface to Xray Data"
    def __init__(self, server='sqlite', dbname=None):
        self.dbname = dbname
        self.server = server
        self.tables = None
        self.engine = None
        self.session = None
        self.conn    = None
        self.metadata = None
        self.pvs = {}
        self.restoring_pvs = []
        if dbname is not None or server=='mysql':
            self.connect(dbname, server=server)

    def connect(self, dbname, server='sqlite'):
        "connect to an existing database"
        if server == 'sqlite':
            if not os.path.exists(dbname):
                raise IOError("Database '%s' not found!" % dbname)

            if not isXrayDB(dbname):
                raise ValueError("'%s' is not an Xray Database file!" % dbname)

        self.dbname = dbname
        self.engine = make_engine(dbname, server)
        self.conn = self.engine.connect()
        self.session = sessionmaker(bind=self.engine)()
        self.metadata =  MetaData(self.engine)
        self.metadata.reflect()
        tables = self.tables = self.metadata.tables

        try:
            clear_mappers()
        except:
            pass
        mapper(ChantlerTable,      tables['Chantler'])
        mapper(WaasmaierTable,     tables['Waasmaier'])

    def close(self):
        "close session"
        self.session.flush()
        self.session.close()

    def query(self, *args, **kws):
        "generic query"
        return self.session.query(*args, **kws)

    def _get_foreign_keyid(self, table, value, name='name',
                           keyid='id', default=None):
        """generalized lookup for foreign key arguments
    table: a valid table class, as mapped by mapper.
    value: can be one of the following
         table instance:  keyid is returned
         string:          'name' attribute (or set which attribute with 'name' arg)
            a valid id
            """
        if isinstance(value, table):
            return getattr(table, keyid)
        else:
            if isinstance(value, (str, unicode)):
                xfilter = getattr(table, name)
            elif isinstance(value, int):
                xfilter = getattr(table, keyid)
            else:
                return default
            try:
                query = self.query(table).filter(
                    xfilter==value)
                return getattr(query.one(), keyid)
            except (IntegrityError, NoResultFound):
                return default

        return default

