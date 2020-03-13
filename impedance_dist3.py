import numpy as np
from numpy.linalg import inv
import argparse, json, yaml, sys

# module globals
rho=100.0   # conductivity of the earth
f=60.0
w = 2*np.pi*f

mu0 = 4*np.pi*np.power(10.0,-7.0)

m_per_mi = 1609.34

# Earth impedance
rg = m_per_mi*mu0*w/8 # convert to units of Ohms/mi
# xg = m_per_mi*mu0*w/(2*np.pi)*np.log(658.5*np.sqrt(rho/f))
xg = 7.93402    

def map_phases(phasing,rem_zeros=False):
    sphasing = sorted(phasing)

    if rem_zeros:
        phases = sphasing
    else:
        phases = 'ABCN'

    mapping = {}

    for i,p in enumerate(phases):
        if p in sphasing:
            j = sphasing.index(p)
            mapping[i] = phasing[j]

    return mapping


def get_phases(phasing):
    return [(p in phasing) for p in 'ABC']


def kron(A):
    B = np.copy(A[0:-1,0:-1])
    n = B.shape[0]    
    
    for r in xrange(0,n):
        for c in xrange(0,n):
            B[r,c] -= A[r,n]*A[n,c]/A[n,n]
    
    return B


def seqz(Z):
    if len(Z) == 1:
        return (Z,Z)

    # Bergen & Vittal p.485
    a = np.exp(1j*2.0*np.pi/3.0)

    Zs0 = (Z[0,0] + Z[1,1] + Z[2,2])/3.0
    # Zsp = (Z[0,0] + a*Z[1,1] + a*a*Z[2,2])/3.0
    # Zsn = (Z[0,0] + a*a*Z[1,1] + a*Z[2,2])/3.0    

    Zm0 = (Z[0,1] + Z[1,2] + Z[0,2])/3.0
    # Zmp = (Z[0,1] + a*Z[1,2] + a*a*Z[0,2])/3.0
    # Zmn = (Z[0,1] + a*a*Z[1,2] + a*Z[0,2])/3.0

    Zp = Zs0 - Zm0
    Z0 = Zs0 + 2*Zm0

    return (Z0,Zp)

def calc_z_prim(D,gmr,rac):       
    n = D.shape[0]

    Z = np.zeros(D.shape,dtype=complex)

    for r in xrange(0,n):
        for c in xrange(0,n):
            if r == c:
                ri = rac[r]
                Z[r,c] = ri + rg + 1j*0.12134*(np.log(1./gmr[r]) + xg)
            else:
                Z[r,c] = rg + 1j*0.12134*(np.log(1./D[r][c]) + xg)

    return Z # return primitive Z matrix


def calc_oh_line(line):       
    # convert GMR to units of m

    # if configs[self.ph_cond].Imax:
    #     Imax_ph = configs[self.ph_cond].Imax
    n = 3

    spacing = line['spacing']
    Zfull = np.zeros((n,n),dtype=complex)

    phase_map = map_phases(line['phasing'],rem_zeros=True)

    gmr = [line['ph_cond']['GMR']]*n
    rac = [line['ph_cond']['r']]*n

    if 'neu_cond' in line:
        n += 1
        gmr.append(line['neu_cond']['GMR'])
        rac.append(line['neu_cond']['r'])
    
    D = np.zeros((n,n))

    for r in xrange(0,n):
        for c in xrange(0,n):
            if (r not in phase_map) or (c not in phase_map):
                continue

            if r == c:
                continue

            pr = phase_map[r]
            pc = phase_map[c]
            
            rp = min(pr,pc)
            cp = max(pr,pc)

            D[r,c] = spacing[rp][cp]

    Zfull = calc_z_prim(D,gmr,rac)
    has_phase = get_phases(line['phasing'])

    if 'neu_cond' in line:
        Z = kron(Zfull) 
        has_phase = has_phase[:-1]
    else:
        Z = Zfull

    # if 'Imax' in configs[line['ph_cond']]:     
    #     line['Imax'] = [Imax_ph]*len(Z)
    
    # line['phases'] = has_phase

    line['Z'] = Z
    line['D'] = D

    (z0,z1) = seqz(line['Z'])
    line['z0'] = z0
    line['z1'] = z1     

    return line


# GMR units: ft
# diameter units: in

def calc_ug_line(line,rem_zeros=False,equiv_neu=False):
    # fill in ph_cond & neu_cond
    cable_cfg = line['cable_cfg']
    # Imax_ph = cable_cfg['Imax']
    core_cfg = line['core_cfg']
    D = line['spacing']

    # calculate neutral GMR & spacings
    # don't need to make any changes to ph_cond_cfg
    D_outer = cable_cfg['D_outer'] # outer diameter (up to outside of neutral strands, not incl. jacket) in in.

    n_ph = len(line['phasing'].replace('N',''))
    n_neu = n_ph

    # need to change GMR of neu_cond_cfg
    if cable_cfg['neu_type'] == 'concentric':
        # import ipdb; ipdb.set_trace()
        strand_cfg = cable_cfg['strand_cfg']
        k = float(cable_cfg['n_strands'])
        gmr_strand = strand_cfg['GMR']
        D_strand = strand_cfg['D'] # diameter of neutral strands in in.
        R_strand = (D_outer - D_strand)/24.0 # distance from cable center to neutral strand centers

        gmr_cable_neu = np.power(gmr_strand * k * np.power(R_strand,k-1.0), 1.0/k)
        # print 'Neutral GMR', gmr_cable_neu
        r_cable_neu = strand_cfg['r']/k
    else: # neutral type is 'Tape'
        # import pdb; pdb.set_trace()

        T = cable_cfg['tape_thick'] # tape thickness in mils
        gmr_cable_neu = ((D_outer/2.0) - (T/2000.0))/12.0
        r_cable_neu = 1.0636e9*cable_cfg['rho_m']/(D_outer*T)


        if 'neutral' in line:
            gmr_neu = configs[line['neutral']]['GMR']
            r_neu = configs[line['neutral']]['r']
            n_neu += 1

    phase_map = map_phases(line['phasing'],rem_zeros=True)
    has_phase = get_phases(line['phasing'])

    Zpp = np.zeros((n_ph,n_ph),dtype=complex)
    Zpn = np.zeros((n_ph,n_neu),dtype=complex)
    Znn = np.zeros((n_neu,n_neu),dtype=complex)
    Dpp = np.zeros((n_ph,n_ph)) # Dpp = Dnn
    Dpn = np.zeros((n_ph,n_neu))
    Dnn = np.zeros((n_neu,n_neu))

    # move this into zprim or calc_line
    # calculate distance matrix & partitioned impedance matrix for phases
    for r in xrange(0,n_ph):
        for c in xrange(r,n_ph):
            if r == c: # self-impedance of phase conductors
                Dpp[r,c] = core_cfg['GMR']
            else: # mutual impedance of phase conductors
                 # Mutual Impedances
                # r <= c
                rp = phase_map[r]
                cp = phase_map[c]

                Dpp[r,c] = Dpp[c,r] = D[rp][cp]

    Zpp = calc_z_prim(Dpp,[core_cfg['GMR']]*n_ph,[core_cfg['r']]*n_ph)

    gmr = [gmr_cable_neu]*n_neu
    rac = [r_cable_neu]*n_neu
    # calculate distance matrix & partitioned impedance matrix for neutrals
    for r in xrange(0,n_neu):
        if r == (n_neu - 1) and 'neutral' in line:
            gmr[r] = gmr_neu
            rac[r] = r_neu

        for c in xrange(r,n_neu):
            if r == c: # self-impedance of phase conductors
                if c == (n_neu - 1) and 'neutral' in line:
                    Dnn[r,c] = gmr_neu
                else:
                    Dnn[r,c] = gmr_cable_neu

            else: # mutual impedance of phase conductors
                 # Mutual Impedances
                # r <= c
                rp = phase_map[r]
                cp = phase_map[c]

                Dnn[r,c] = Dnn[c,r] = D[rp][cp]
           
    Znn = calc_z_prim(Dnn,gmr,rac)

    # calculate distance matrix & partitioned impedance matrix for neutrals
    for r in xrange(0,n_ph):
        for c in xrange(r,n_neu):
            if r == c: # GMD of cable phase to its own concentric neutral
                if cable_cfg['neu_type'].lower() == 'concentric':
                    Dpn[r,c] = R_strand   
                else:
                    Dpn[r,c] = gmr_cable_neu

                Zpn[r,c] = rg + 1j*0.12134*(np.log(1./Dpn[r,c]) + xg)   
                # print 'Cab Neu GMR', Dpn[r,c]
                # print 'xi', np.log(1./Dpn[r,c]), 'xg', xg
                # print 'Xpn', Zpn[r,c].imag.round(6)
            elif c == (n_neu - 1): # GMD of cable phase to external concentric neutral
                rp = phase_map[r]
                cp = phase_map[c]

                Dpn[r,c] = D[rp][cp]
                Zpn[r,c] = rg + 1j*0.12134*(np.log(1./Dpn[r,c]) + xg)
            else: # GMD of cable phase to another cable's concentric neutral
                 # Mutual Impedances
                # r <= c
                rp = phase_map[r]
                cp = phase_map[c]

                if cable_cfg['neu_type'].lower() == 'concentric':
                    if equiv_neu:
                        Dpn[r,c] = np.sqrt(D[rp][cp]*D[rp][cp] + R_strand*R_strand)
                    else:
                        Dpn[r,c] = np.power(np.power(D[rp][cp],k) - np.power(R_strand,k), 1/k)
                else:
                    Dpn[r,c] = np.sqrt(D[rp][cp]*D[rp][cp] + gmr_cab_neu*gmr_cab_neu)
             
                Zpn[r,c] = rg + 1j*0.12134*(np.log(1./Dpn[r,c]) + xg)

            if c < n_ph:
                Zpn[c,r] = Zpn[r,c]

    # do Kron reduction
    line['Z'] = Zpp - np.dot(np.dot(Zpn,inv(Znn)),np.transpose(Zpn)) # Zpn = Znp
    line['D'] = M = np.bmat( [[Dpp, Dpn], [np.transpose(Dpn), Dnn]] )

    # import ipdb; ipdb.set_trace()
    (z0,z1) = seqz(line['Z'])
    line['z0'] = z0
    line['z1'] = z1     
    
    # if 'Imax' in cable_cfg:
    #     line['Imax'] = [Imax_ph]*n_ph
    
    line['phases'] = has_phase
    # line['rg'] = round(rg,6)
    # line['xg'] = round(xg,6)
    return line

def calc_xf(xf):
    xf['Z'] = np.zeros((3,3),dtype=complex)

    # calculate the base impedance
    zb = 1e3*xf['Vpri']*xf['Vpri']/xf['Smax']

    for k in xrange(0,3):
        xf['Z'][k,k] = zb*(xf['Rpu'] + 1j*xf['Xpu']) 

    return xf