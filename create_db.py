'''
This script creates an SQLite3 database of fundamental X-ray fluorescence
parameters as compiled by W.T. Elam, B.D. Ravel and J.R. Sieber, published in
Radiation Physics and Chemistry, 63 (2), 121 (2002). The database is published
by NIST at http://www.cstl.nist.gov/acd/839.01/xrfdownload.html
'''

import io
import json
import os
import sqlite3
import sys
from string import maketrans


    
def add_f0Waasmaier(dest, append=True):
    """add f0 data from Waasmaier and Kirfel"""
    source = 'f0_WaasKirf.dat'

    if os.path.exists(dest) and not append:
        raise IOError('File "%s" already exists -- cannot add f0 data')

    conn = sqlite3.connect(dest)
    c = conn.cursor()
    c.execute(
        '''create table f0Waasmaier (id integer,
        atomic_number integer, element text, ion text,
        offset real, scale text, exponents text)
        ''')
    
    f = open(source)
    lines = f.readlines()
    if 'Elastic Photon-Atom Scatt' not in lines[1]:
        raise RuntimeError('Source file not recognized for f0_WaasKirf data')

    strip_ion = maketrans('0123456789+-', ' '*12)
    id = 0
    while lines:
        line = lines.pop(0)        
        if line.startswith('#S '):
            id += 1
            #print [s for s in line[3:].split()]
            zstr, ion = [s.strip() for s in line[3:].split()]
            atno = int(zstr)
            for i in range(3):
                line = lines.pop(0)
            words = [float(w.strip()) for w in line.split()]
            off   = words[5]
            scale = json.dumps(words[:5])
            expon = json.dumps(words[6:])

            elem = ion.translate(strip_ion).strip()
            c.execute('insert into f0Waasmaier values (?,?,?,?,?,?,?)',
                      (id, atno, elem, ion, off, scale, expon))

    conn.commit()
    c.close()


def add_chantler(dest, append=True):
    """add f' / f'' data from Chantler"""

    if os.path.exists(dest) and not append:
        raise IOError('File "%s" already exists -- cannot add f0 data')

    conn = sqlite3.connect(dest)
    c = conn.cursor()
    c.execute(
        '''create table chantler (id integer,
        element text, sigma_mu real, mue_f2 real, density real,
        energy text, f1 text, f2 text, mu_photo text,
        mu_incoh text, mu_total text)
        ''')

    dirname = 'chantler'
    args = '(%s)' % ','.join(['?']*11)
    
    nelem = 92
    for z in range(1, nelem+1):
        fname = os.path.join(dirname, '%2.2i.dat' % z)
        lines = open(fname, 'r').readlines()

        # line 1: take symbol and density only
        words = lines[0][1:-1].split()
        words.pop()
        density = float(words.pop())
        elem = words[0].replace(':','')
        
        # line 2: take sigma_mu
        words = lines[1][1:-1].split()
        sigma_mu = float(words.pop())
        
        # line 3: take mue_f2
        words = lines[2][1:-1].split()
        mue_f2 = float(words.pop())

        en, f1, f2, mu_photo, mu_incoh, mu_total = [], [], [], [], [], []
        for line in lines:
            if line.startswith('#'):
                continue
            words = [float(w) for w in line[:-1].split()]
            en.append(words[0])
            f1.append(words[1]  - z)
            f2.append(words[2])
            mu_photo.append(words[3])
            mu_incoh.append(words[4])
            mu_total.append(words[5])

        c.execute('insert into chantler values %s' % args,
                  (z, elem, sigma_mu, mue_f2, density,
                   json.dumps(en), json.dumps(f1), json.dumps(f2), 
                   json.dumps(mu_photo), json.dumps(mu_incoh), 
                   json.dumps(mu_total)))
        
    conn.commit()
    c.close()

def create_database(source, dest, overwrite=False, silent=False):
    if not os.path.isfile(source):
        if silent:
            return
        raise IOError('File "%s" does not exist' % source)
    if os.path.isfile(dest) and overwrite:
        os.remove(dest)
    if os.path.exists(dest):
        if silent:
            return
        raise IOError('File "%s" already exists. Use "-f" to overwrite' % dest)

    with io.open(source, encoding='ascii') as f:
        lines = f.readlines()
        if 'Elam, Ravel, Sieber' not in lines[0]:
            raise RuntimeError('Source file not recognized')
        while lines[0].startswith('/'):
            lines.pop(0)

    conn = sqlite3.connect(dest)
    c = conn.cursor()

    c.execute(
        '''create table elements (atomic_number integer, element text,
        molar_mass real, density real)
        '''
        )
    current_edge_id = 0
    c.execute(
        '''create table xray_levels (id integer, element text, iupac_symbol
        text, absorption_edge real, fluorescence_yield real, jump_ratio real)
        '''
        )
    current_line_id = 0
    c.execute(
        '''create table xray_transitions (id integer, element text,
        iupac_symbol text, siegbahn_symbol text, initial_level text,
        final_level text, emission_energy real, intensity real)
        '''
        )
    current_ck_id = 0
    c.execute(
        '''create table Coster_Kronig
        (id integer, element text, initial_level text, final_level text,
        transition_probability real, total_transition_probability real)
        '''
        )
    current_photo_id = 0
    c.execute(
        '''create table photoabsorption (id integer, element text,
        log_energy text, log_photoabsorption text,
        log_photoabsorption_spline text)
        '''
        )
    current_scatter_id = 0
    c.execute(
        '''create table scattering (id integer, element text, log_energy text,
        log_coherent_scatter text, log_coherent_scatter_spline text,
        log_incoherent_scatter text, log_incoherent_scatter_spline text)
        '''
        )

    while lines:
        line = lines.pop(0)
        if line.startswith('Element'):
            sym, num, mw, rho = line.split()[1:]
            c.execute(
                'insert into elements values (?,?,?,?)', (num, sym, mw, rho)
                )
            current_element = sym
        elif line.startswith('Edge'):
            current_edge_id += 1
            label, energy, yield_, jump = line.split()[1:]
            el = current_element
            c.execute(
                'insert into xray_levels values (?,?,?,?,?,?)',
                (current_edge_id, el, label, energy, yield_, jump)
                )
            current_edge = label
        elif line.startswith('  Lines'):
            while True:
                if lines[0].startswith('    '):
                    current_line_id += 1
                    line = lines.pop(0)
                    iupac, siegbahn, energy, intensity = line.split()
                    start, end = iupac.split('-')
                    el = current_element
                    c.execute(
                        'insert into xray_transitions values (?,?,?,?,?,?,?,?)',
                        (current_line_id, el, iupac, siegbahn, start, end,
                        energy, intensity)
                        )
                else:
                    break
        elif line.startswith('  CK '):
            temp = line.split()[1:]
            ck = []
            while temp:
                (i,j), temp = temp[:2], temp[2:]
                ck.append((i,j))
            if lines[0].startswith('  CKtotal'):
                temp = lines.pop(0).split()[1:]
                ck_total = []
                while temp:
                    (i,j), temp = temp[:2], temp[2:]
                    ck_total.append((i,j))
            else:
                ck_total = ck
            for i, j in zip(ck, ck_total):
                current_ck_id += 1
                (so, p), tp = i[:], j[1]
                c.execute(
                    '''insert into Coster_Kronig
                    values (?,?,?,?,?,?)''',
                    (current_ck_id, current_element, current_edge, so, p, tp)
                    )
        elif line.startswith('Photo'):
            current_photo_id += 1
            energy = []
            photo = []
            spline = []
            while lines[0].startswith('    '):
                temp = [float(i) for i in lines.pop(0).split()]
                energy.append(temp[0])
                photo.append(temp[1])
                spline.append(temp[2])
            c.execute(
                'insert into photoabsorption values (?,?,?,?,?)',
                (current_photo_id, current_element, json.dumps(energy),
                json.dumps(photo), json.dumps(spline))
                )
        elif line.startswith('Scatter'):
            current_scatter_id += 1
            energy = []
            cs = []
            css = []
            ics = []
            icss = []
            while lines[0].startswith('    '):
                temp = [float(i) for i in lines.pop(0).split()]
                energy.append(temp[0])
                cs.append(temp[1])
                css.append(temp[2])
                ics.append(temp[3])
                icss.append(temp[4])
            c.execute(
                'insert into scattering values (?,?,?,?,?,?,?)',
                (current_scatter_id, current_element, json.dumps(energy),
                json.dumps(cs), json.dumps(css), json.dumps(ics),
                json.dumps(icss))
                )

    conn.commit()

    c.close()


if __name__ == '__main__':
    try:
        import argparse
    except ImportError:
        raise RuntimeError(
            'argparse module not found.\n'
            'Install argparse or update to python-2.7 or >=python-3.2.'
            )
    parser = argparse.ArgumentParser(
        description='export the Elam ascii file to an SQLite database "dest"'
        )
    parser.add_argument('source')
    parser.add_argument('dest')
    parser.add_argument('-f', '--force', action='store_true')
    parser.add_argument('-s', '--silent', action='store_true')
    args = parser.parse_args()

    create_database(
        args.source, args.dest, overwrite=args.force, silent=args.silent
        )

    add_f0Waasmaier(args.dest, append=True)
    add_chantler(args.dest, append=True)
