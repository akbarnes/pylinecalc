import toml, logging, carsons
import numpy as np
from copy import deepcopy
from pprint import pprint

# logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)


# self impedance from Kersting 4.41 (zii)
calc_zs = lambda r, GMR: r + 0.0953 + 0.12134j*(7.93402 - np.log(GMR))

# mutual impedance from Kersting 4.42 (zij)
calc_zm = lambda Dij: 0.0953 + 0.12134j*(7.93402 - np.log(Dij))

# Kersting 5.9
calc_ps = lambda Sii, RDi: 11.17689*np.log(Sii/RDi)

# Kersting 5.10
calc_pm = lambda Sij, Dij: 11.17689*np.log(Sij/Dij)

capacitance = lambda Pij: np.linalg.inv(Pij)


def shunt_admittance(Cij, f=60):
    return np.pi*f*Cij


def calc_zprim(D, gmr, rac):       
    n = D.shape[0]
    Z = np.zeros(D.shape, dtype=complex)

    for r in range(0, n):
        for c in range(0, n):
            if r == c:
                Z[r,c] = calc_zs(rac, gmr)
            else:
                Z[r,c] = calc_zm(D[r,c])

    return Z # return primitive Z matrix


# kron reduction
def kron (Zij, Zin, Znn): 
    if Zin is None or Znn is None:
        return Zij

    if np.isscalar(Znn):
        return Zij - np.outer(Zin,Zin)/Znn # single neutral form
    
    return Zij - np.linalg.solve(np.matmul(Zin, Zin.T), Znn)


# for interactive use, use short_keys=True
# for calling in a gui use short_keys=False
def load_configs(filename, long_line_keys=False, duplicate_line_keys=True):
    raw_configs = toml.load(filename)
    configs = {}

    for table_name, table in raw_configs.items():
        # logging.debug(f'Updating {table_name}')
        configs[table_name] = {}

        for record_key,record in table.items():
            if table_name =='lines':
                if long_line_keys and 'name' in record:
                    configs[table_name][record['name']] = deepcopy(record)

                    if duplicate_line_keys:
                        configs[table_name][record_key] = deepcopy(record)
                else:
                    configs[table_name][record_key] = deepcopy(record)
            else:
                configs[table_name][record_key] = deepcopy(record)

                if 'name' in record:
                    configs[table_name][record['name']] = deepcopy(record)


    # logging.debug('Setting conductors & spacings')
    conductors = configs['conductors']
    spacings = configs['spacings']
    cables = configs['cables']
    
    for cfg_name,cfg in configs['cables'].items():
        cfg['core'] = conductors[cfg['core']]

        if 'strand' in cfg:
            cfg['strand'] = conductors[cfg['strand']]        

    for cfg_name,cfg in configs['lines'].items():

        # for OH lines 
        if 'phase' in cfg and isinstance(cfg['phase'], str):
            cfg['phase'] = conductors[cfg['phase']]

        if 'neutral' in cfg and isinstance(cfg['neutral'], str):
            cfg['neutral'] = conductors[cfg['neutral']]

        # for UG lines
        if 'cable' in cfg and isinstance(cfg['cable'], str):
            cfg['cable'] = cables[cfg['cable']]

        if 'spacings' in cfg and isinstance(cfg['spacings'], str):
            cfg['spacings'] = spacings[cfg['spacings']]

     
            
    return configs

def overhead_line(line): 
    phasing = 'abcn'
    
    if phasing in line:
        phasing = line['phasing'].lower()
    
    phases = phasing.replace('n', '')
    
    phase = line['phase']
    D = line['spacings']
  
    Zij = np.zeros((3,3), dtype=complex)

    for i,pi in enumerate('abc'):
        for j,pj in enumerate('abc'):
            pij = f'{pi}{pj}'

            if i == j and pi in phases:
                Zij[i,i] = calc_zs(phase['r'], phase['GMR'])          
            elif pij in D:
                Zij[i,j] = Zij[j,i] = calc_zm(D[pij])
            elif pi in D and pj in D:
                xi,yi = D[pi]
                xj,yj = D[pj]
                
                Dij = np.sqrt((xi - xj)**2 + (yi - yj)**2)
                Zij[i,j] = Zij[j,i] = calc_zm(Dij)

    if 'neutral' not in line:
        return Zij, None, None

    Zin = np.zeros((3,1), dtype=complex)

    for i,pi in enumerate('abc'):
        pin = f'{pi}n'

        if pin in D:
            Zin[i] = calc_zm(D[pin])

    neutral = line['neutral']    
    znn = calc_zs(neutral['r'], neutral['GMR'])  
    return Zij, Zin, znn

class LineModel:
    def __init__(self, conductors):
        self.resistance = {}
        self.geometric_mean_radius = {}
        self.wire_positions = {}
        self.phases = {}

        for phase, (r, gmr, (x, y)) in conductors.items():
            self.resistance[phase] = r
            self.geometric_mean_radius[phase] = gmr
            self.wire_positions[phase] = (x, y)
            self.phases = sorted(list(conductors.keys()))

def carsons_overhead_line(line, f=60): 
    # build the conductors
    phases = line['phasing'].lower()
    conds = {}

    for p in phases:
        if p in 'abc':
            r = line['phase']['r']
            gmr = line['phase']['GMR']
        else:
            r = line['neutral']['r']
            gmr = line['neutral']['GMR']  

        x,y = line['spacings'][p]        
        conds[p.upper() ] = r, gmr, (x,y)

    lm = LineModel(conds)
    pprint(lm.__dict__)
    return carsons.convert_geometric_model(lm) 

def overhead_line_shunt(line):
    phasing = 'abcn'
    
    if phasing in line:
        phasing = line['phasing'].lower()
        
    phasing = line['phasing'].lower()
    phases = phasing.replace('n', '')
    
    phase = line['phase']
    D = line['spacings']
  
    Pij = np.zeros((3,3), dtype=complex)
       
    for i,pi in enumerate('abc'):
        for j,pj in enumerate('abc'):
            if i == j and pi in phases:
                Sii = 2*D[pi][1]
                Pij[i,i] = calc_ps(Sii, phase['D']/12)
            elif pi in D and pj in D:
                xi, yi = D[pi]
                xj, yj = D[pj]
                
                Dij = np.sqrt((xi - xj)**2 + (yi - yj)**2)
                Sij = np.sqrt((xi - xj)**2 + (yi + yj)**2)
                Pij[i,j] = calc_pm(Sij, Dij)

    if 'neutral' not in line:
        return Pij, None, None

    Pin = np.zeros((3,1), dtype=complex)
    xn, yn = D['n']

    for i,pi in enumerate('abc'):
        if pi in D:
            xi, xj = D[pi]
            
            Din = np.sqrt((xi - xn)**2 + (yi - yn)**2)
            Sin = np.sqrt((xi - xn)**2 + (yi + yn)**2)
            Pin[i,0] = calc_pm(Sin, Din)

    neutral = line['neutral']    
    Snn = 2*yn
    pnn = calc_ps(Snn, neutral['D']/12)  
    return Pij, Pin, pnn


def underground_line(line):
    cable = line['cable']
    has_concentric_neutral = False
    
    if cable['neutral_type'].lower() == 'concentric':
        has_concentric_neutral = True
    
    pprint(line)
    phasing = 'abcn'
    
    if 'phasing' in line:
        phasing = line['phasing'].lower()
        
    phases = phasing.replace('n','')
    print(f'phasing = {phasing}, phases = {phases}')
    D = line['spacings']

    if 'ca' in D:
        D['ac'] = D['ca'] # make life a bit easier for creating the spacing matrix

    D_outer = cable['D_outer']
    n_ph = n_neu = len(phases)
    
    if has_concentric_neutral:
        strand = cable['strand']
        k = cable['n_strands']
        gmr_strand = strand['GMR']
        D_strand = strand['D'] # diameter of neutral strands in in.
        R_strand = (D_outer - D_strand)/24 # distance from cable center to neutral strand centers
        
        gmr_cable_neu = np.power(gmr_strand * k * np.power(R_strand, k-1), 1/k)
        r_cable_neu = strand['r']/k
    else: # neutral type is tape
        T = cable['tape_thickness']
        gmr_cable_neu = (D_outer/2 - T/2000)/12
        r_cable_neu = 1.0636e9*cable['rho_m']/(D_outer*T)
        
    Zpn = np.zeros((n_ph,n_neu), dtype=complex)
    Dpp = np.zeros((n_ph,n_ph)) # Dpp = Dnn
    Dpn = np.zeros((n_ph,n_neu))
    Dnn = np.zeros((n_neu,n_neu))

    core = cable['core']
    logging.debug(core)

    # move this into zprim or calc_line
    # calculate distance matrix & partitioned impedance matrix for phases
    for r in range(0, n_ph):
        for c in range(r, n_ph):
            if r == c: # self-impedance of phase conductors
                Dpp[r,c] = core['GMR']
            else: # mutual impedance of phase conductors
                rp = phases[r]
                cp = phases[c]

                Dpp[r,c] = Dpp[c,r] = D[f'{rp}{cp}']

    Zpp = calc_zprim(Dpp, core['GMR'], core['r'])

    # calculate distance matrix & partitioned impedance matrix for neutrals
    for r in range(0, n_neu):
        for c in range(r, n_neu):
            if r == c: # self-impedance of cable neutral
                Dnn[r,c] = gmr_cable_neu
            else: # mutual impedance of cable neutral
                rp = phases[r]
                cp = phases[c]

                Dnn[r,c] = Dnn[c,r] = D[f'{rp}{cp}']
           
    Znn = calc_zprim(Dnn, gmr_cable_neu, r_cable_neu)
    
    # calculate distance matrix for phase-neutrals
    for r in range(0, n_ph):
        for c in range(r, n_neu):
            if r == c: # GMD of cable phase to its own concentric neutral
                if has_concentric_neutral:
                    Dpn[r,c] = R_strand
                else:
                    Dpn[r,c] = gmr_cable_neu
            else: # GMD of cable phase to another cable's concentric neutral
                 # Mutual Impedances
                # r <= c
                rp = phases[r]
                cp = phases[c]
                prc = f'{rp}{cp}'

                if has_concentric_neutral:
                    Dpn[r,c] = Dpn[c,r] = np.sqrt(D[prc]**2 + R_strand**2)
                else:
                    Dpn[r,c] = Dpn[c,r] = np.sqrt(D[prc]**2 + gmr_cable_neu**2)
             
    # calculate partitioned impedance matrix for phase-neutrals
    for r in range(0, n_ph):
        for c in range(r, n_neu):
            if r == c: # GMD of cable phase to its own concentric neutral
                Zpn[r,c] = calc_zm(Dpn[r,c])  
            else: # GMD of cable phase to another cable's concentric neutral
                 # Mutual Impedances
                # r <= c
                rp = phases[r]
                cp = phases[c]
                prc = f'{rp}{cp}'
                Zpn[r,c] = Zpn[c,r] = calc_zm(Dpn[r,c])

    return Zpp, Zpn, Znn, Dpp, Dpn, Dnn


# this assumes XLPE dielectric
def underground_line_shunt(line):
    cable = line['cable']
    has_concentric_neutral = False 
    core = cable['core']

    
    if cable['neutral_type'].lower() == 'concentric':
        has_concentric_neutral = True
    
    phasing = 'abcn'
    
    if phasing in line:
        phasing = line['phasing'].lower()
        
    phases = phasing.replace('n','')
    
    Yij = np.zeros((3,3), dtype=complex)

    for i,pi in enumerate('abc'):
        if pi in phases:
            RDi = core['D']/24 # radius of the center conductor in ft
            logging.debug(f'Radius of central conductor: RDi = {RDi:0.6f} ft')                

            D_outer = cable['D_outer']

            if has_concentric_neutral:
                k = cable['n_strands']
    
                strand = cable['strand']

                # diameter of neutral strands in in.
                # this is D_strand in underground_line()
                D_strand = strand['D']
                RDs = strand['D']/24
            
                # distance from cable center to neutral strand centers
                # this is R_strand in underground_line()
                Rb = (D_outer - strand['D'])/24 
                
                logging.debug(f'Radius of circle passing through strands: Rb = {Rb:0.6f} ft')
                logging.debug(f'Radius of strands: RDi = {RDs:0.6f} ft')                                
            
                Yij[i,i] = 77.3619j/(np.log(Rb/RDi) - np.log(k*RDs/Rb)/k) # Siemens/mile
            else:
                T = cable['tape_thickness']                
                Rb = (D_outer/2 - T/2000)/12
                logging.debug(f'Radius of circle passing through tape: Rb = {Rb:0.6f} ft')
                
                Yij[i,i] = 77.3619j/np.log(Rb/RDi) # Siemens/mile
                
    return Yij


calc_zself = lambda Z: (Z[0,0] + Z[1,1] + Z[2,2])/3
calc_zmut = lambda Z: (Z[0,1] + Z[1,2] + Z[0,2])/3

calc_z0 = lambda zs, zm: zs + 2*zm
calc_z1 = lambda zs, zm: zs - zm

def sequences(Z):
    zs = calc_zself(Z)
    zm = calc_zmut(Z)
    
    z0 = calc_z0(zs, zm)
    z1 = calc_z1(zs, zm)
    
    return z0, z1

    