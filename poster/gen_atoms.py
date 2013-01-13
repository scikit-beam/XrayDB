from xraydb import xrayDB
db = xrayDB()

def extract(dat, key):
    val = dat.get(key, [' '])[0]
    if isinstance(val, float) and val > 1:
        val = int(round(val))
    return str(val)

elnames = ('', 'hydrogen' , 'helium', 'lithium' , 'beryllium', 'boron' ,
           'carbon' ,'nitrogen' , 'oxygen', 'fluorine' , 'neon' ,'sodium',
           'magnesium', 'aluminum' , 'silicon', 'phosphorus' , 'sulfur',
           'chlorine' , 'argon', 'potassium' , 'calcium', 'scandium' ,
           'titanium','vanadium' , 'chromium','manganese' , 'iron',
           'cobalt' , 'nickel','copper' , 'zinc','gallium' ,
           'germanium','arsenic' , 'selenium','bromine' , 'krypton',
           'rubidium' , 'strontium','yttrium' , 'zirconium','niobium' ,
           'molybdenum','technetium' , 'ruthenium','rhodium' ,
           'palladium','silver' , 'cadmium','indium' , 'tin', 'antimony' ,
           'tellurium','iodine' , 'xenon','cesium' , 'barium','lanthanum' ,
           'cerium','praseodymium', 'neodymium', 'promethium' ,
           'samarium','europium' , 'gadolinium', 'terbium' ,
           'dysprosium','holmium' , 'erbium','thulium' ,
           'ytterbium','lutetium' , 'hafnium','tantalum' , 'tungsten',
           'rhenium' , 'osmium','iridium' , 'platinum','gold' ,
           'mercury','thallium' , 'lead','bismuth' , 'polonium', 'astatine',
           'radon','francium' , 'radium','actinium' ,
           'thorium','protactinium', 'uranium','neptunium' , 'plutonium',
           'americium', 'curium', 'berkelium', 'californium')

table=r"""\newcommand{\%(texsym)s}{{%%
\begin{minipage}{67mm}

\vspace{1mm}

{\Huge{\hspace{1mm} {\textbf{%(sym)s}} \hfill \hfil{\textbf{%(iz)s}} \hspace{1mm}}} %%

\vspace{6mm}

{\Huge{\hfill {\Name{%(name)s}} \hfill}}

\vspace{6mm}

{\Large{
\begin{tabular*}{67mm}%%
{@{\hspace{5pt}}{r}@{\extracolsep{\fill}}r@{\extracolsep{\fill}}r}%%
{\BRed{%(k)s}}  & {\bf{%(ka1)s}} &  {\bf{%(kb1)s}} \\%%
{\BBlue{%(l1)s}} & %(lb3)s & %(lb4)s \\%%
{\BBlue{%(l2)s}} & %(lb1)s & %(lg1)s \\%%
{\BRed{%(l3)s}} & {\bf{%(la1)s}} & {\bf{%(lb2)s}} \\%%
{\BBlue{%(m5)s}} & %(ma)s & %(mb)s \\%%
\multicolumn{3}{@{\hspace{1pt}}c}{ }\\%%
\end{tabular*}
\vfill}}
\end{minipage}}}
"""

for iz in range(1, 99):
    sym = db.symbol(iz)

    dat = {'iz': iz, 'sym': sym, 'name': elnames[iz],
           'texsym': "Elem%s" % sym}

    edges = db.xray_edges(iz)
    lines = db.xray_lines(iz)
    dat['k']   = extract(edges, 'K')
    dat['ka1'] = extract(lines, 'Ka1')
    dat['kb1'] = extract(lines, 'Kb1')

    dat['l1']  = extract(edges, 'L1')
    dat['lb3'] = extract(lines, 'Lb3')
    dat['lb4'] = extract(lines, 'Lb4')

    dat['l2']  = extract(edges, 'L2')
    dat['lb1'] = extract(lines, 'Lb1')
    dat['lg1'] = extract(lines, 'Lg1')

    dat['l3']  = extract(edges, 'L3')
    dat['la1'] = extract(lines, 'La1')
    dat['lb2'] = extract(lines, 'Lb2,15')

    dat['m5'] = extract(edges, 'M5')
    dat['ma'] = extract(lines, 'Ma')
    dat['mb'] = extract(lines, 'Mb')

    # print dat.keys()
    print '%% ', dat['name']
    print table % dat

highz  = ((99, 'Es', 'einsteinium'),
         (100, 'Fm', 'fermium'),
         (101, 'Md', 'mendelevium'),
         (102, 'No', 'nobelium'),
         (103, 'Lr', 'lawrencium'))

for iz, sym, name in highz:
    dat = {'iz': iz, 'sym': sym, 'name': name,
           'texsym': "Elem%s" % sym,
           'k': '', 'ka1': '', 'kb1': '', 'l1' : '', 'lb3': '', 'lb4': '',
           'l2' : '', 'lb1': '', 'lg1': '', 'l3' : '', 'la1': '', 'lb2':
           '', 'm5': '', 'ma': '', 'mb': ''}

    print '%% ', dat['name']
    print table % dat

print '%% Key'
print table % {'iz': 'Z', 'sym': 'Symbol', 'name': 'name',
               'texsym': "ElemKey",
               'k':   r'$\mathbf{K}$ edge',
               'ka1': r'$\mathbf{K_{\alpha_1}}$',
               'kb1': r'$\mathbf{K_{\beta_1}}$',
               'l1':  r'$\mathrm{L_{\rm I}}$ edge',
               'lb3': r'$\mathrm{L_{\beta_3}}$',
               'lb4': r'$\mathrm{L_{\beta_4}}$',
               'l2':  r'$\mathrm{L_{\rm II}}$ edge',
               'lb1': r'$\mathrm{L_{\beta_1}}$',
               'lg1': r'$\mathrm{L_{\gamma_1}}$',
               'l3':  r'$\mathbf{L_{\rm III}}$ edge',
               'la1': r'$\mathbf{L_{\alpha_1}}$',
               'lb2': r'$\mathbf{L_{\beta_2}}$',
               'm5':  r'$\mathrm{M_{\rm V}}$ edge',
               'ma':  r'$\mathrm{M_{\alpha}}$',
               'mb':  r'$\mathrm{M_{\beta}}$'}

