import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import glob
import scipy.spatial as spatial

# Debugging. 
os.chdir('/mnt/c/Users/User/Documents/Academia/Fall_2023/Research/Simulations/WSL/CO_v3_esf/tests/dwfn')

# region: Variables. 
ry2ev = 13.60569
bohr2a = 0.529177249
step_size = 0.01/bohr2a
wfn = np.array([])
elph = np.array([])
dwfn_wfn = np.array([])
dwfn_elph = np.array([])
# endregion 

# region: Functions. 
def read_wfn(file):
    global alat
    global nat
    global avec
    global atyp
    global vb_max 
    global eigs

    with h5py.File(file, 'r') as wfn:        
        alat = wfn['/mf_header/crystal/alat'][()]   # Bohr. 
        nat = wfn['/mf_header/crystal/nat'][()]
        avec = np.einsum('ij->ji', wfn['mf_header/crystal/avec'][:])     # Lattice vectors, unit of alat. (lattice vector, axis).
        atyp = wfn['mf_header/crystal/atyp'][:]     # Atomic number of each atom. 
        vb_max = wfn['mf_header/kpoints/ifmax'][:][0, 0] - 1 # Assume all values in this array are the same for non-metals. (Since I'm only interested in them atm).
        wfn_out = wfn['wfns/coeffs'][:]
        eigs = wfn['mf_header/kpoints/el'][:][0, 0, :]      # Size should be (500,) or what ever is the number of bands used for getting WFN.h5. 
    
    # Cha
    
    return np.squeeze(np.vectorize(complex)(wfn_out[..., 0], wfn_out[..., 1]), axis=1)

def get_kgrid():        # Get kgrid. 
    global kgrid 

    kgrid = np.loadtxt('./kgrid.out', skiprows=2)

    # Additional checking to get the right size. 
    if (len(kgrid.shape) == 1):
        kgrid = kgrid.reshape((1, kgrid.shape[0]))

    kgrid = kgrid[:, [0, 1, 2]]

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
    
    elph_folded = np.zeros((nmodes, nbnd_ph, nkpts_ph, nbnd_ph, nkpts_ph), dtype='c16')

    # From (nmodes, nbnd, nq, nbnd, nk) -> (nmodes, nbnd, nk, nbnd, nkp). 
    # Strategy: Find k+q -> k_r. This will replace the 3rd dimension. 

    for mode in range(nmodes):
        for bnd in range(nbnd_ph):
            for qpt in range(nqpts_ph):
                for bndp in range(nbnd_ph):
                    for kptp in range(nkpts_ph):

                         elph_folded[mode, bnd, get_k_plus_q_idx(kptp, qpt), bndp, kptp] = elph[mode, bnd, qpt, bndp, kptp]

def read_elph():
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

def get_from_wfn():
    '''
    Will get the dwfn from wfn. 
    '''
    global dwfn_wfn

    # Allocate. 
    dwfn_wfn = np.zeros((nbnd_ph, nbnd_ph), dtype='c16')

    # Assemble. 
    for bnd in range(nbnd_ph):
        for bndp in range(nbnd_ph):

            diff = wfn_posdel[bndp, :] - wfn_pos[bndp, :]
            deriv = diff/step_size
            bra = np.conjugate(wfn_pos[bnd, :])
            dwfn_wfn[bnd, bndp] = np.dot(bra, deriv)

def get_from_elph():
    global dwfn_elph

    dwfn_elph = elph_folded[0, :, 0, :, 0]            # (mode, i, k, j, kp). We choose the first mode (x direction for C). 

    # Divide by the energy term. 
    for bnd in range(nbnd_ph):
        for bndp in range(nbnd_ph):

            ebnd = eigs[bnd]
            ebndp = eigs[bndp]

            if np.abs(ebnd - ebndp) < 1e-4:
                dwfn_elph[bnd, bndp] = 0.0
            else:
                dwfn_elph[bnd, bndp] /= (ebndp - ebnd)      # Second hellman feynman theorem. 

# endregion

# region: Main. 

# Read data. 
get_kgrid()
wfn_pos = read_wfn('WFN_1.15.h5')
wfn_posdel = read_wfn('WFN_1.16.h5')
read_elph()

# Assemble data. 
get_from_wfn()
get_from_elph()

x = 1

# Print and compare. Maybe plot. 
fig = plt.figure(figsize=(12, 6))

ax = fig.add_subplot(1, 2, 1)
im = ax.pcolormesh(np.abs(dwfn_wfn), cmap='RdBu')
fig.colorbar(im, ax=ax)
ax.set_title('Calculated using WFN')

ax = fig.add_subplot(1, 2, 2)
im = ax.pcolormesh(np.abs(dwfn_elph), cmap='RdBu')
fig.colorbar(im, ax=ax)
ax.set_title('Calculated using elph')

plt.show()

# endregion. 