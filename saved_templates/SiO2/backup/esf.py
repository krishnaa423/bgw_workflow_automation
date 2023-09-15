# Other general modules. 
import numpy as np
import matplotlib.pyplot as plt 
from ase import Atoms
from ase.io import espresso
import h5py 
import os
import time 
import re
import scipy.linalg as linalg
from scipy import spatial
import pyparsing as pp
import glob as glob

# Modules for pygsl. 
from pygsl import errno
from pygsl import errors
from pygsl import math
from pygsl import multiminimize
from pygsl import vector


# Will strive to keeps units at (A, eV) for (length, energy). Phonon frequency in THz. 
# Assumptions. q-grid is the same as k-grid. For now, no interpolation scheme. For testing. Also no symmetries are used, so no irz is assumed.
os.chdir('/mnt/c/Users/User/Documents/Academia/Fall_2023/Research/Simulations/WSL/CO_v3_esf/8-esf')

# Set plot style. 
plt.style.use('bmh')



# region: Variables.

# region: optimization. 
step_size, max_tol = 0.01, 1e-4
max_steps = 20
max_force = 0.02
niter = max_steps       # Just an alias for max_steps. 
iter = 0                # Current iteration index. Zero based.
status = 0
dim = -1
params = 0
optim_obj = 0
my_eps = 1e-5
# endregion

# region: STE parameters.
xct_state_index = 0         # 0 based. 
# endregion

# region: Units. 
ry2ev = 13.60569
bohr2a = 0.529177249
# endregion

# region: WFN.
alat = 0.0                  # Lattice parameter in bohr. 
nat  = 2                  # In case of CO. 
avec = np.array([])         # (lattice vector, axis). 
ntyp = np.array([])         # Array of atomic numbers. 
xctd_struct = np.array([])    # Shape: (max_steps+1, nat, 3). Units: A. 
nkpts_dft = 1
vb_max  =  4             # valence band max 0 indexed from bottom.                   
# endregion

# region: dft_force.
dft_force = np.array([])  # Shape: (nat, 3). Units: eV/A.  
# endregion

# region: elph.
nbnd_elph_min = 0         # 0 indexed. 
nbnd_elph_max = 19        # 0 indexed. 
nbnd_ph = 20             # Might have to manually set this. Or automate from workflow. 
nkpts_ph = 1
nqpts_ph = 1
nmodes = 6
elph = np.array([])    # Shape: (nbnd_ph, nbnd_php, nk, nmodes, nq). Units: eV/A. 
elph_folded = np.array([])   # Shape: (nmodes, nbnd_ph, nk, nbnd_php, nkp). Units: eV/A.
# endregion

# region: kgrid.
kgrid = np.array([])   # Size: (nk, 3). Units: None. 
# endregion

# region: gw_eig.
nbnd_gw_min = 0         # 0 indexed. 
nbnd_gw_max = 9        # 0 indexed. 
nbnd_gw = 10            # Might have to manually set this. Or automate from workflow. 
nkpts_gw = 1
gw_eig = np.array([])  # Shape: (nbnd_gw, nk). Units: eV. 
# endregion

# region: xct_vec.
ns = 1
nbnd_cond = 5
nbnd_val = 5
nkpts_abs = 1
ncvks = nbnd_cond*nbnd_val*nkpts_abs
xct_vec = np.array([])   # Shape: (nbnd_cond, nbnd_val, nk): Units: None. 
# endregion

# region: xct_e
xct_e = 0.0              # Single val based on xct_state_index. 
# endregion

# region: hbse.
hbse = np.array([])    # Shape: (nbnd_cond, nbnd_val, nk, nbnd_condp, nbnd_valp, nkp). Units: eV. 
K = np.array([])       # Shape: (nbnd_cond, nbnd_val, nk, nbnd_condp, nbnd_valp, nkp). Units: eV. 
# endregion

# region: ESF.
# Calculations.
P = np.array([])        # Shape: (nmodes, nbnd_ph, nk, nbnd_php, nkp). Units: A-1.
dK = np.array([])       # Shape: (nmodes, nbnd_cond, nbnd_val, nk, nbnd_condp, nbnd_valp, nkp). Units: eV/A.
dGW = np.array([])      # Shape: (nmodes, nbnd_cond, nbnd_val, nk, nbnd_condp, nbnd_valp, nkp). Units: eV/A.
dhbse = np.array([])    # Shape: (nmodes, nbnd_cond, nbnd_val, nk, nbnd_condp, nbnd_valp, nkp). Units: eV/A.
dOmega = np.array([])   # Shape: (nmodes). Units: eV/A.  nmodes = (nat, 3), although might be (3, nat) if you think in Fortan order.  

# Storing results.
xctd_e = np.array([])           # Shape: (niter). Units: eV.
xctd_force = np.array([])       # Shape: (niter, nmodes). Units: eV/A. 
xctd_struct = np.array([])      # Shape: (niter, nmodes). Units: A. nmodes = (nat, 3).flatten() in C order.
# endregion

# endregion.

# region: Functions.

# region: Helper functions. 
def get_double_flip(x, axes=[1, 3]):
    return np.flip(
        np.flip(
            x,
            axis=axes[1]
        ),
        axis=axes[0]
    )

def get_gw_cond_bnd(bnd_idx):   # bnd_idx is zero indexed. 
    idx =  (vb_max - nbnd_gw_min) + 1 + bnd_idx
    return idx

def get_gw_val_bnd(bnd_idx):
    idx = (vb_max - nbnd_gw_min) - bnd_idx
    return idx

def get_elph_cond_bnd(bnd_idx):    
    idx = (vb_max-nbnd_elph_min) + 1 + bnd_idx
    return idx

def get_elph_val_bnd(bnd_idx):     
    idx = (vb_max - nbnd_elph_min) - bnd_idx
    return idx

def get_gw_from_elph_bnd(elph_bnd_idx): 
    idx = (elph_bnd_idx + nbnd_elph_min) - nbnd_gw_min
    return idx

def get_k_plus_q_idx(k, q): # Note, the kpts are not in the order you might think. 
    
    kdtree = spatial.KDTree(kgrid)

    sum_vec = kgrid[k, :] + kgrid[q, :]

    # Bring vector back to Brillouin zone. 
    for idx, value in enumerate(sum_vec):
        if value >= 1.0: sum_vec[idx] -= 1.0
        if value < 0.0: sum_vec[idx] += 1.0

    # Ask what k-point index does this correspond to. 
    return kdtree.query(sum_vec)[1]

def get_elph_folded():
    global elph_folded

    if iter == 0:
        elph_folded = np.zeros((nmodes, nbnd_ph, nkpts_ph, nbnd_ph, nkpts_ph), dtype='c16')

    # From (nmodes, nbnd, nq, nbnd, nk) -> (nmodes, nbnd, nk, nbnd, nkp). 
    # Strategy: Find k+q -> k_r. This will replace the 3rd dimension. 

    for mode in range(nmodes):
        for bnd in range(nbnd_ph):
            for qpt in range(nqpts_ph):
                for bndp in range(nbnd_ph):
                    for kptp in range(nkpts_ph):

                         elph_folded[mode, bnd, get_k_plus_q_idx(kptp, qpt), bndp, kptp] = elph[mode, bnd, qpt, bndp, kptp]

def get_K_from_gw_eig():
    global K

    # Shape of hbse: (nbnd_cond, nbnd_val, nkpts_abs, nbnd_cond, nbnd_val, nkpts_abs). Shape of gw: (nbnd_gw, nkpts_gw)
    # Shape of K is the same as the shape of hbse. 

    K = np.zeros((nbnd_cond, nbnd_val, nkpts_abs, nbnd_cond, nbnd_val, nkpts_abs), dtype='c16')

    # Update K. 
    for bnd_cond in range(nbnd_cond):
        for bnd_val in range(nbnd_val):
            for kpt in range(nkpts_abs):

                # Set the diagonal. Which is just the gap energies. 
                K[bnd_cond, bnd_val, kpt, bnd_cond, bnd_val, kpt] -= gw_eig[get_gw_cond_bnd(bnd_cond), kpt] - gw_eig[get_gw_val_bnd(bnd_val), kpt]

    K += hbse

def is_energy_close(a, b):
    if np.abs(a - b) < 1e-5: 
        return True 
    else: 
        False
# endregion. 


# region: Get inputs.

def get_WFN():          # Get apos -> xctd_struct, nat, nkpts, vb_max, nmodes from WFN.h5
    global alat
    global nat
    global avec
    global atyp
    global xctd_struct
    global vb_max 

    with h5py.File('WFN.h5', 'r') as wfn:
        
        alat = wfn['/mf_header/crystal/alat'][()]   # Bohr. 
        nat = wfn['/mf_header/crystal/nat'][()]
        avec = np.einsum('ij->ji', wfn['mf_header/crystal/avec'][:])     # Lattice vectors, unit of alat. (lattice vector, axis).
        atyp = wfn['mf_header/crystal/atyp'][:]     # Atomic number of each atom. 
        
        # Allocate if it is the first iteration.
        if iter == 0:
            xctd_struct = np.zeros((max_steps+1, nat*3))
            xctd_struct[iter, :] = (wfn['mf_header/crystal/apos'][:]).flatten(order='C')*alat*bohr2a   # Cartesial coordinates of each atom in units of angstrom. (nat, 3)
        vb_max = wfn['mf_header/kpoints/ifmax'][:][0, 0] - 1 # Assume all values in this array are the same for non-metals. (Since I'm only interested in them atm).

def get_dft_force():    # Get df_force from pwscf.out
    global dft_force

    if iter == 0:
        dft_force = np.zeros((nat*3))

    pattern = ...+ pp.Literal('Forces acting on atoms') + ...+ pp.Literal(':\n\n') + ...+ pp.Literal('The non-local')

    parsed_txt = pattern.parse_file('pwscf.out').as_list()

    force_lines = parsed_txt[4].splitlines()

    for atm_index, line in enumerate(force_lines):
        dft_force[atm_index*3:(atm_index+1)*3] =  np.array(list(map(lambda x: float(x), line.split()[-3:])))*ry2ev/bohr2a   # In eV/A units. 

def get_atm_masses():   # Get atomic masses from pwscf.in. atm_masses = (nat). 
    global atm_masses

    # Init size. 
    if iter == 0:
        atm_masses = np.zeros(nat)

    # Set pattern. 
    pattern = ...+ pp.Literal('ATOMIC_SPECIES\n') + pp.Regex(r'.*', re.DOTALL)

    parsed_txt = pattern.parse_file('./pwscf.in').as_list()

    species_lines = pattern.parse_file('./pwscf.in').as_list()[2].splitlines()

    # Add to atm_masses data. 
    for species_idx, species_line in enumerate(species_lines):
        if species_idx < nat:
            atm_masses[species_idx] = species_line.split()[1]
        else:       # Only take into account <nat> lines. 
            break

def get_kgrid():        # Get kgrid. 
    global kgrid 

    kgrid = np.loadtxt('./kgrid.out', skiprows=2)

    # Additional checking to get the right size. 
    if (len(kgrid.shape) == 1):
        kgrid = kgrid.reshape((1, kgrid.shape[0]))

    kgrid = kgrid[:, [0, 1, 2]]

def get_elph():         # Get elph from elph.dat. 
    global nbnd_ph
    global nkpts_ph
    global nmodes
    global nqpts_ph
    global elph

    # Almost always complex flavor for calculations. 
    flavor = 2

    # TODO: Need to read elph band window info from epw.in.

    # Read from elph_*.h5 files. 
    elph_files = glob.glob('./elph_*.h5')
    for kpt_idx, elph_file in enumerate(elph_files):
        with h5py.File(elph_file, 'r') as elph_h5:
            elph_shape = elph_h5['/elph_cart_real'][:].shape

            if kpt_idx == 0:  # Need only do this once. 
                # Shape is (q, nu, k, j, i). 
                nqpts_ph = elph_shape[0]
                nmodes_ph = elph_shape[1]
                nmodes = nmodes_ph
                nkpts_ph = len(elph_files)
                nbnd_ph = elph_shape[3]

                # Init the size of elph. 
                elph = np.zeros((nqpts_ph, nmodes_ph, nkpts_ph, nbnd_ph, nbnd_ph), dtype='c16')

            # Read the elph and assign the matrix. 
            # Create complex array. Assuming the units are in Ry/bohr initially. Converted to eV/A. 
            elph[:, :, kpt_idx, :, :] = np.squeeze(elph_h5['/elph_cart_real'][:] + 1j*elph_h5['/elph_cart_imag'][:], axis=2)*ry2ev/bohr2a

    # Reorder the array. 
    elph = np.einsum('qmkji->miqjk', elph)  # Reorder.    New shape is (nmodes, nbnd_ph, nqpts_ph, nbnd_ph, nkpts_ph). 

    # Calculate the folded matrix. 
    get_elph_folded()       # Changes the global elph_folded variable. Size: (nmodes, nbnd_ph, nkpts_ph, nbnd_ph, nkpts_ph)

def get_gw_eig():       # Get gw_eig from eqp.dat.
    global gw_eig
    global nbnd_gw_min
    global nbnd_gw_max
    global nbnd_gw
    global nkpts_gw

    data = np.loadtxt('./eqp.dat')   # 4 columns (qx, qy, qz, nbnd_gw) for first line, then nbnd_gw lines have (_, bnd_num, gw_eig, _)

    # Get nbnd_gw. 
    nbnd_gw = int(data[0, 3])

    # Get nbnd_gw_min
    nbnd_gw_min = int(data[1, 1]) - 1

    # Get nbdn_gw_max. 
    nbnd_gw_max = int(data[nbnd_gw, 1]) - 1

    # Get nkpts_gw.
    num_rows = data.shape[0]
    nkpts_gw = int(num_rows/(nbnd_gw+1))

    if iter == 0:
        gw_eig = np.zeros((nbnd_gw, nkpts_gw))

    # Extract the gw_energies for each k-point. 
    for kpt_gw in range(nkpts_gw):
        gw_eig[:, kpt_gw] = data[kpt_gw*(nbnd_gw+1) + 1:(kpt_gw+1)*(nbnd_gw+1), 2]

def get_xct_vec_and_eig():      # Get xct_vec, ncvks from eigenvectors.h5. 
    global xct_vec
    global ns       # Have to do this for spin case. 
    global nbnd_cond
    global nbnd_val
    global nkpts_abs
    global xctd_e

    # Init for first iter. 
    if iter == 0:
        xctd_e = np.zeros(max_steps+1)

    xct = h5py.File('eigenvectors.h5', 'r')

    # Only working under TDA approximation. 
    ns = xct['exciton_header/params/ns'][()]
    nbnd_cond = xct['exciton_header/params/nc'][()]
    nbnd_val = xct['exciton_header/params/nv'][()]
    nkpts_abs = xct['exciton_header/kpoints/nk'][()]
    nQ = xct['exciton_header/kpoints/nQ'][()] # We are always going to have nQ = 1. 

    xctd_e[iter] = xct['exciton_data/eigenvalues'][:][xct_state_index]

    # Shape before squeeze: (nk, nc, nv, ns (squeezed), flavor). After squeeze: (nk, nc, nv, flavor). 
    xct_vec = np.squeeze(xct['exciton_data/eigenvectors'][:][0, xct_state_index, ...], axis=3)    
    xct_vec = np.vectorize(complex)(xct_vec[...,0], xct_vec[...,1]) # Create complex array. (nk, nc, nv). 

    # Use einsum to permute to the form in the equations. 
    print(f'Size of xct_vec before einsum: {xct_vec.shape}')
    xct_vec = np.einsum('kcv->cvk', xct_vec)    
    print(f'Size of xct_vec after einsum: {xct_vec.shape}')

def get_hbse_and_K():         # Get hbse, K from hbse_a.dat.
    global hbse
    global xct_vec

    # Read from hbse.h5 file. 
    with h5py.File('hbse.h5', 'r') as hbse_h5:
        
        # TDA approximation and reshape. 
        hbse = np.reshape(
            hbse_h5['/hbse_a'][:]*ry2ev, 
            (nkpts_abs, nbnd_cond, nbnd_val, nkpts_abs, nbnd_cond, nbnd_val, 2),   # New shape: (k, c, v, k', c', v', flavor). 
            order='C'
        )

        # Convert to complex.
        hbse = np.vectorize(complex)(hbse[..., 0], hbse[..., 1])

        # Permute using einsum to convenient order in equations. 
        # hbse = np.einsu('kcvKCV->cvkCVK', hbse) 
        print(f'Size of hbse before einsum: {hbse.shape}')
        hbse = np.einsum('kcvKCV->cvkCVK', hbse) 
        print(f'Size of hbse after einsum: {hbse.shape}')


    # Overwrite the xct_vec for debugging. 
    ncvk = nbnd_cond*nbnd_val*nkpts_abs

    # Calc eigs and evecs. 
    hbse_reshaped = np.reshape(hbse, (ncvk, ncvk), order='C')
    eigs, evecs = linalg.eig(hbse_reshaped)

    eigs_sort = np.argsort(eigs)

    calc_evec = evecs[:, eigs_sort[xct_state_index]]

    xct_vec = calc_evec.reshape((nbnd_cond, nbnd_val, nkpts_abs))

    # Get the K matrix.
    get_K_from_gw_eig()

# endregion.

# region: Calc intermediates and output. 

def calc_P():
    global P

    P = np.zeros(
        (
        nmodes,
        nbnd_ph,
        nkpts_ph,
        nbnd_ph,
        nkpts_ph
        ),
        dtype='c16'
    )
    
    # TODO: Can be vectorized. 
    # P is (nmodes, nbnd_ph, nk, nbnd_php, nkp). elph_folded is (nmodes, nbnd_ph, nk, nbnd_php, nkp). 
    # gw_eig is (nbnd_gw, nk). 
    for mode in range(nmodes):
        for bnd in range(nbnd_ph):
            for kpt in range(nkpts_gw):
                for bndp in range(nbnd_ph):
                    for kptp in range(nkpts_gw):
                        
                        # Calculate index in gw array based on index in elph array. 
                        gw_bnd_idx = get_gw_from_elph_bnd(bnd)
                        gw_bndp_idx = get_gw_from_elph_bnd(bndp)


                        # Skip if calculated band index is not within gw band range. 
                        if gw_bnd_idx < 0 or gw_bnd_idx > (nbnd_gw - 1):
                            continue
                        if gw_bndp_idx < 0 or gw_bndp_idx > (nbnd_gw - 1):
                            continue


                        if not is_energy_close(gw_eig[gw_bnd_idx, kpt], gw_eig[gw_bndp_idx, kptp]):
                            try:
                                energy_diff = gw_eig[gw_bndp_idx, kptp] - gw_eig[gw_bnd_idx, kpt]
                            except Exception as e:
                                print(f'Exception occured creating energy diff for P. Info: {e}')

                            P[mode, bnd, kpt, bndp, kptp] = elph_folded[mode, bnd, kpt, bndp, kptp]/energy_diff

def calc_dK():
    global dK

     # Init with zeros. 
    dK = np.zeros(
        (
        nmodes,
        nbnd_cond,
        nbnd_val,
        nkpts_abs,
        nbnd_cond,
        nbnd_val,
        nkpts_abs
        ),
        dtype='c16'
    )

    # Calculate some row limits of P.
    P_cond_start = get_elph_cond_bnd(0)
    P_cond_end = get_elph_cond_bnd(nbnd_cond - 1)
    P_val_end = get_elph_val_bnd(0)
    P_val_start = get_elph_val_bnd(nbnd_val - 1)

    # Calc values using einsum. 
    # 4 terms. In the order listed in the (2003, Beigi) ESF paper. \

    # 1
    dK += np.einsum(
        'Mjqck,jvqCVK->McvkCVK',
        np.conjugate(
            P[
                :, 
                P_cond_start:P_cond_end+1,
                :,
                P_cond_start:P_cond_end+1,
                :
            ]
        ),
        K
    )

    # 2
    dK += np.einsum(
        'Mjqvk,cjqCVK->McvkCVK',
        get_double_flip(P[
            :, 
            P_val_start:P_val_end+1,
            :,
            P_val_start:P_val_end+1,
            :
        ]),
        K
    )

    # 3
    dK += np.einsum(
        'MjqCK,cvkjVq->McvkCVK',
            P[
            :, 
            P_cond_start:P_cond_end+1,
            :,
            P_cond_start:P_cond_end+1,
            :
        ],
        K
    )

    # 4
    dK += np.einsum(
        'MjqVK,cvkCjq->McvkCVK',
        get_double_flip(np.conjugate(
                P[
                :, 
                P_val_start:P_val_end+1,
                :,
                P_val_start:P_val_end+1,
                :
            ],
        )),
        K
    )
 
def calc_dGW():
    global dGW

    # Init with zeros. 
    dGW = np.zeros(
        (
        nmodes,
        nbnd_cond,
        nbnd_val,
        nkpts_abs,
        nbnd_cond,
        nbnd_val,
        nkpts_abs
        ),
        dtype='c16'
    )

    # TODO: Can be vectorized. 
    # elph_folded has shape (nmodes, nbnd_ph, nk, nbnd_ph, nk). Set the diagonal entries. 
    for mode in range(nmodes):
        for bnd_cond in range(nbnd_cond):
            for bnd_val in range(nbnd_val):
                for kpt in range(nkpts_abs):

                    # Set the matrix entries.  
                    dGW[
                        mode, 
                        bnd_cond, 
                        bnd_val,
                        kpt,
                        bnd_cond, 
                        bnd_val,
                        kpt
                    ] = elph_folded[mode, get_elph_cond_bnd(bnd_cond), kpt, get_elph_cond_bnd(bnd_cond), kpt] 
                    - elph_folded[mode, get_elph_val_bnd(bnd_val), kpt, get_elph_val_bnd(bnd_val), kpt]

def calc_dhbse():
    global dhbse

    dhbse = dGW + dK

def calc_dOmega():
    global xctd_force

    if iter == 0:
        xctd_force = np.zeros((max_steps+1, nmodes), dtype='c16')

    # Calculate the mat/vec size. 
    ncvk = nbnd_cond*nbnd_val*nkpts_abs

    # Reshape. 
    xct_vec_reshaped = xct_vec.reshape((ncvk), order='C')
    dhbse_reshaped = dhbse.reshape((nmodes, ncvk, ncvk), order='C')

    # Einsum for product. A'*dhbse*A. 
    dOmega = np.einsum(
        'i,kij,j->k',
        np.conjugate(xct_vec_reshaped),
        dhbse_reshaped,
        xct_vec_reshaped
    )

    gwbse_force = -dOmega 

    # Calc the gw contrib to force for debugging. 
    dhgw = dGW 
    dhgw_reshaped = dhgw.reshape((nmodes, ncvk, ncvk), order='C')

    # Einsum for product. A'*dhbse*A. 
    dOmega_gw = np.einsum(
        'i,kij,j->k',
        np.conjugate(xct_vec_reshaped),
        dhgw_reshaped,
        xct_vec_reshaped
    )

    gw_force = -dOmega_gw 



    # Set the force for the current iteration.
    print(f'DFT force: {np.real(dft_force)}') 
    print(f'GW force alone: {np.real(gw_force)}') 
    print(f'BSE force force: {np.real(gwbse_force - gw_force)}') 
    print(f'GW-BSE force alone: {np.real(gwbse_force)}') 
    print(f'Total force: {np.real(gwbse_force + dft_force)}') 
    xctd_force[iter, :] = gwbse_force + dft_force

def calc_xctd_force():
    
    calc_P()
    calc_dK()
    calc_dGW()
    calc_dhbse()
    calc_dOmega()

def calc_new_xctd_struct():
    global xctd_struct

    # Set the structure for the next iteration. 
    xctd_struct[iter+1, :] = optim_obj.getx()     # Get the current location. 

# endregion.

# region: Functions for ESD. 

def is_within_my_eps(x, y):

    # Check if any one is out of error bound. 
    for idx, val in enumerate(x):
        if np.abs(val, y[idx]) > my_eps:
            return False
    
    return True     # Only gets here if the above fails. 

def return_e(v, params):
    
    # Run calculation if necessary. Then return exciton energy. 
    if not is_within_my_eps(xctd_struct[iter, :], v):
        get_e_and_esf()
    else:
        return xctd_e[iter, :]

def return_esf(v, params):
    
    # Run calculation if necessary. Then return ESF. 
    if not is_within_my_eps(xctd_struct[iter, :], v):
        get_e_and_esf()
    else:
        return xctd_force[iter, :]

def return_e_and_esf(v, params):
    # Run calculation if necessary. Then return ESF. 
    if not is_within_my_eps(xctd_struct[iter, :], v):
        get_e_and_esf()
    else:
        return xctd_e[iter, :], xctd_force[iter, :]

def setup_optim():
    global optim_obj        # Modify the global object. 

    dim = nmodes
    my_fdf = multiminimize.gsl_multimin_function_fdf(return_e, return_esf, return_e_and_esf, params, dim)  # Set the functions. 
    optim_obj = multiminimize.vector_bfgs(my_fdf, dim) # Set the minimizer state. 

# region: WORKFLOW

# region: Timer variables. 
sched_type = 'WSL'   # WSL or slurm. 
is_time_logging = True   # Want to record how much time each one takes? 
# endregion

# region: Timer class. 
class Logger:

    def __init__(self):
        self.start = 0.0
        self.stop = 0.0

    def time_start(self):
        self.start = time.time()

    def time_stop(self):
        self.stop = time.time()

    def time_elapsed(self):
        return self.stop - self.start
local_logger = Logger()
global_logger = Logger()
# endregion

# region: Workflow timer functions. 
def print_flush(msg):
    print(msg, flush=True)

def run_cmd(cmd):
    global sched_type
    global local_logger
    global global_logger

    # Starting a command. 
    print_flush(f'Starting : {cmd}')
    exit_code = -1

    # Run the command. 
    if is_time_logging: local_logger.time_start()

    if sched_type == 'WSL':
        exit_code = os.system(f'{cmd}')
    elif sched_type == 'slurm':
        exit_code = os.system(f'sbatch --wait {cmd}')

    if is_time_logging: local_logger.time_stop()

    # See if it finished properly. 
    if exit_code == 0:  # No error. 
        # Done with a command. 
        if is_time_logging:
            print_flush(f'Done with: {cmd}. Took {local_logger.time_elapsed()} seconds.\n\n')
        else:
            print_flush(f'Done with: {cmd}.\n\n')
    else:
        # Error. 
        if is_time_logging:
            print_flush(f'Error finishing: {cmd}. Exited with code {exit_code}. Took time {local_logger.time_elapsed()}.\n\n')
            global_logger.time_stop()
            print_flush(f'Total time elapsed for the workflow run is: {global_logger.time_elapsed()} seconds.') 

        else:
            print_flush(f'Error finishing: {cmd}. Exited with code {exit_code}.\n\n')
            global_logger.time_stop()
            print_flush(f'Total time elapsed for the workflow run is: {global_logger.time_elapsed()} seconds.') 

        os._exit(exit_code)     # Exit with the error code. 
# endregion

# region: Start global timer
global_logger.time_start()
# endregion

def get_e_and_esf():        # With global structure info. 
    
    # Run the whole workflow. 

    # Write current struct to input files: 
    # 1--- 1.1-scf/pwscf.in
    # 2--- 1.2-ph/pwbands.in. Also 7---1.4-dos/pwbands.in.
    # 3--- 2.1-wfn/pwscf.in
    # 4--- 2.2-wfnq/pwbands.in
    # 5--- 3.1-wfn_fi/pwbands.in
    # 6--- 3.2-wfnq_fi. 
    # Data from xctd_struct[iter, :] In Fortran order (3, nat) for the positions in angstrom. 

    if iter != 0:       # Need not update cell parameters for first iteration. 
        # Files. 
        files = [
            '../1.1-scf/pwscf.in',
            '../1.2-ph/pwbands.in',
            '../2.1-wfn/pwbands.in',
            '../2.2-wfnq/pwbands.in',
            '../3.1-wfn_fi/pwbands.in',
            '../3.2-wfnq_fi/pwbands.in'
        ]

        # Create substitution text. 
        subst_txt = 'CELL_POSITIONS angstrom\n'
        apos = xctd_struct[iter, :].reshape(nat, 3, order='C')
        for row in apos: subst_txt + f'{row[0]} {row[1]} {row[2]}\n'

        # TODO: grammar. Also hook subst text action.
        pattern = ... + pp.Literal('CELL_POSITIONS') + ... + pp.Char(pp.srange('[A-Z]')) + pp.Regex(r'.*', re.DOTALL)

        # Replace in files. 
        for file in files:
            parsed_txt = pattern.parse_file(file).as_list()

            # Create transformed_text.
            transformed_txt = ''.join([
                parsed_txt[0],
                subst_txt,
                parsed_txt[3],
                parsed_txt[4]
            ]) 

            # Write to file. 
            with open(file, 'w') as f: f.write(transformed_txt)

    # Run all code. 
    # region: 1.1-scf
    os.chdir('../1.1-scf')
    run_cmd('./jobscript_pwscf.run')
    os.chdir('../')
    # endregion

    # region: 1.2-ph
    os.chdir('./1.2-ph')
    run_cmd('./jobscript_pwbands.run')
    run_cmd('./jobscript_ph.run')
    # run_cmd('./jobscript_q2r.run')
    # run_cmd('./jobscript_matdyn.run')
    os.chdir('../')
    # endregion

    # region: 1.3-epw
    os.chdir('./1.3-epw')
    run_cmd('./jobscript_epw.run')
    os.chdir('../')
    # endregion

    # region: 2.1-wfn
    os.chdir('./2.1-wfn')
    run_cmd('./jobscript_kgrid.run')
    run_cmd('./jobscript_pwbands.run')
    run_cmd('./jobscript_pw2bgw.run')
    # run_cmd('./jobscript_parabands.run')
    os.chdir('../')
    # endregion

    # region: 2.2-wfnq
    # os.chdir('./2.2-wfnq')
    # run_cmd('./jobscript_kgrid.run')
    # run_cmd('./jobscript_pwbands.run')
    # run_cmd('./jobscript_pw2bgw.run')
    # run_cmd('./jobscript_parabands.run')
    # os.chdir('../')
    # endregion


    # region: 3.1-wfn_fi
    # os.chdir('./3.1-wfn_fi')
    # run_cmd('./jobscript_kgrid.run')
    # run_cmd('./jobscript_pwbands.run')
    # run_cmd('./jobscript_pw2bgw.run')
    # run_cmd('./jobscript_parabands.run')
    # os.chdir('../')
    # endregion


    # region: 3.2-wfnq_fi
    # os.chdir('./3.2-wfnq_fi')
    # run_cmd('./jobscript_kgrid.run')
    # run_cmd('./jobscript_pwbands.run')
    # run_cmd('./jobscript_pw2bgw.run')
    # os.chdir('../')
    # endregion

    # region: 4-epsilon
    os.chdir('./4-epsilon')
    run_cmd('./jobscript_epsilon.run')
    os.chdir('../')
    # endregion


    # region: 5.1-sigma
    os.chdir('./5.1-sigma')
    run_cmd('./jobscript_sigma.run')
    os.chdir('../')
    # endregion

    # region: 6-kernel
    os.chdir('./6-kernel')
    run_cmd('./jobscript_kernel.run')
    os.chdir('../')
    # endregion


    # region: 7.1-absorption
    os.chdir('./7.1-absorption')
    run_cmd('./jobscript_absorption.run')
    # run_cmd('./jobscript_plotxct.run')          # If you want to plot the exciton too. 
    os.chdir('../')
    # endregion


    # region: 7.2-assemble-kernel
    os.chdir('./7.2-assemble-kernel')
    run_cmd('./jobscript_assemble_kernel.run')
    os.chdir('../')
    # endregion

    # Finally switch to current directory. 
    os.chdir('./8.esf')

    # Load/Update current process stack heap/variables from the above run. 
    get_WFN()
    get_dft_force()
    get_elph()
    get_gw_eig()
    get_hbse_and_K()
    get_xct_vec_and_eig()

    # Calculate esf data. 
    calc_xctd_force()

# endregion

# endregion.

# endregion. 

# region: Debugging.

# region: Test reading
get_WFN()
get_dft_force()
get_atm_masses()
get_kgrid()
get_elph()
get_gw_eig()
get_xct_vec_and_eig()
get_hbse_and_K()
# endregion

# region: Test hbse and elph. 
def test_hbse():    # Calculate the eigenvalues and eigenvectors. 
    global hbse

    ncvk = nbnd_cond*nbnd_val*nkpts_abs

    # Calc eigs and evecs. 
    hbse_reshaped = np.reshape(hbse, (ncvk, ncvk), order='C')
    eigs, evecs = linalg.eig(hbse_reshaped)

    eigs_sort = np.argsort(eigs)

    calc_evec = evecs[:, eigs_sort[xct_state_index]]
    read_evec = np.reshape(xct_vec, (ncvk), order='C')

    norm_calc_evec = np.sum(np.abs(calc_evec)**2)
    norm_read_evec = np.sum(np.abs(read_evec)**2)
    


    # Print values to compare. 
    print(f'Calc HBSE eigenvalues: \n{eigs}', flush=True)
    print(f'Calc HBSE eigenvalues sorted: \n{np.sort(np.real(eigs))}', flush=True)
    print(f'Calc HBSE eigenvalues sort indices: \n{eigs_sort}', flush=True)
    print(f'Read HBSE eigenvalue: \n{xctd_e[iter]}', flush=True)
    print(f'Calc HBSE eigenvector: \n{calc_evec}', flush=True)
    print(f'Calc HBSE eigenvector norm: {norm_calc_evec}', flush=True)
    print(f'Read HBSE eigenvector: \n{read_evec}', flush=True)
    print(f'Read HBSE eigenvector norm: {norm_read_evec}', flush=True)


    # Test read vector eigenvector property. 
    debug_read_eig = np.sqrt(np.sum(np.abs(np.matmul(hbse_reshaped, read_evec))**2))
    debug_calc_eig = np.sqrt(np.sum(np.abs(np.matmul(hbse_reshaped, calc_evec))**2))

    print(f'\n\nLength of  hbse*read_evec: {debug_read_eig}')
    print(f'Length of  hbse*calc_evec: {debug_calc_eig}')


def test_elph():

    # Cant think of a test for now. Just plotting it to see patterns.

    # elph_folded shape: (nmodes, nbnd_ph, nkpts_ph, nbnd_ph, nkpts_ph)  

    elph_folded_reshaped = elph_folded.reshape((nmodes, nbnd_ph*nkpts_ph, nbnd_ph*nkpts_ph))

    fig = plt.figure()
    ax = fig.add_subplot(2, 2, 1)
    im = ax.imshow(np.real(elph_folded_reshaped[0, :, :]))
    fig.colorbar(im, ax=ax)
    ax = fig.add_subplot(2, 2, 2)
    im = ax.imshow(np.real(elph_folded_reshaped[1, :, :]))
    fig.colorbar(im, ax=ax)
    ax = fig.add_subplot(2, 2, 3)
    im = ax.imshow(np.real(elph_folded_reshaped[2, :, :]))
    fig.colorbar(im, ax=ax)
    ax = fig.add_subplot(2, 2, 4)
    im = ax.imshow(np.real(elph_folded_reshaped[3, :, :]))
    fig.colorbar(im, ax=ax)
    plt.show()


# test_hbse()  
# test_elph()   
# endregion.


# region: Get force only. Run full if needed.  
# get_e_and_esf()

# Get force. 
calc_xctd_force()

print(f'Force is: {np.real(xctd_force[iter, :])}.', flush=True)
print(f'Energy is: {xctd_e[iter]}.', flush=True)
# endregion.

# Skip exit if testing above. 
exit(0)

# endregion.


# region: Main loop

# region: Init
setup_optim()
# endregion

# region: Iterate.
while True:
    # Iterate.
    status = optim_obj.iterate()

    # Check was step was okay. 
    if status != errno.GSL_SUCCESS:
        break

    # Calculate new positions if it was okay.
    calc_new_xctd_struct()

    # Get convergence criteria. 
    # print(f'Gradient at iteration {iter}: {s.gradient()}.')
    status = multiminimize.test_gradient(optim_obj.gradient(), max_tol)

    # See if criteria was reached. 
    if (status == errno.GSL_SUCCESS):
        print(f'Minimum force at {xctd_force[iter, :]} after {iter} iterations.')
        print(f'Minimum energy is: {xctd_e[iter, :]}.')
        break
    else:
        print(f'Current minimum energy estimate for iteration {iter}: {xctd_e[iter, :]}.')

    iter += 1
    # Quit if needed. 
    if not ((status == errno.GSL_CONTINUE) and (iter < max_steps)):
        break
# endregion

# endregion.