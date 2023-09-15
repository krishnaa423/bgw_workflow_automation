import struct 
import os
from ase import Atoms
from ase.io import read, write, xsf
import numpy as np
import matplotlib.pyplot as plt
import h5py 

# For debugging
os.chdir('/mnt/c/Users/User/Documents/Academia/Summer_2023/Research/Simulations/WSL/CO_v3_esf/test_elph')

# region: Create CO struct. 
nat = 2
struct_formula = 'CO'
cell = [7, 7, 7, 90, 90, 90] # a, b, c, \alpha, \beta, \gamma.
positions = [        # Units in Angstroms. 
    (0.0, 0.0, 0.0),
    (1.15, 0.0, 0.0),           # Bond length is 1.15 Angstrom. 
]
struc = Atoms(
    struct_formula,
    cell = cell,
    positions= positions
)
struc.center()         # move the molecule to the center. 
# endregion

# region: Read dvscf perturbation grid. 

# fft grid. 
wfn_fft_grid = np.array([36, 36, 36])
dvscf_fft_grid = np.array([72, 72, 72])
size_of_dtype = 16 # 8 for double, and 2 for complex number. 

# Size of files. 
# dvscf file is direct access with 3*nat records (all the perturbations corresponding to a q point).
dvscf_record_len_bytes = int(np.prod(dvscf_fft_grid)*size_of_dtype)
dvscf_record_len_double = int(dvscf_record_len_bytes/8)
dvscf_record_len_complex = int(dvscf_record_len_double/2)
npert = nat*3
dvscf = np.zeros((npert, dvscf_record_len_complex), dtype='c16')

# Read all pertubations for the first (and only) q point.
with open('./CO.dvscf_q1', 'rb') as f:
    for pert in range(npert):
        buffer = f.read(dvscf_record_len_bytes)
        dvscf_row = np.frombuffer(buffer, dtype='f8', count=dvscf_record_len_double)
        dvscf_row = dvscf_row.reshape((dvscf_record_len_complex, 2))
        dvscf_row = np.vectorize(complex)(dvscf_row[..., 0], dvscf_row[..., 1])    # It is a complex number now. 
        dvscf[pert, :] = dvscf_row

dvscf = dvscf.reshape((npert, dvscf_fft_grid[0], dvscf_fft_grid[1], dvscf_fft_grid[2]))
dvscf = np.einsum('mzyx->mxyz', dvscf)   # Reorder from fortran to c order.

# Save xsf file. 
# plot_pert = 1         # In case ony interested in one. 

# Create perturbation files for each one. 
for pert in range(npert):
    with open(f'co_{pert+1}.xsf', 'w') as f:
        xsf.write_xsf(f, [struc], np.abs(dvscf[pert, :]))

print('Done writing xsf files')

# endregion


# Visualize and test. 
# -- Seems good! Based on figures from xcrysden. 

# region: Do the same for wfn data.
bnd_idx = 1

with h5py.File('WFN.h5', 'r') as f:

    wfns = f['/wfns/coeffs'][:]          # Read only the band that we need. 
    wfns = np.vectorize(complex)(wfns[..., 0], wfns[..., 1])
    wfns = wfns.squeeze()

    gvecs = f['/wfns/gvecs'][:]                     # Read the gvectors. 
    # Assume all the cell info is already obtained from before.         

    celvol = f['/mf_header/crystal/celvol'][()]

nbnds = wfns.shape[0]
total_gvecs = int(gvecs.shape[0])

# Bring to FFT grid zone. 
gvecs_folded = np.zeros_like(gvecs)
gvecs_folded[:, 0] = gvecs[:, 0] % wfn_fft_grid[0]
gvecs_folded[:, 1] = gvecs[:, 1] % wfn_fft_grid[1]
gvecs_folded[:, 2] = gvecs[:, 2] % wfn_fft_grid[2]

# Space for FFT. 
wfn_ffts = np.zeros((nbnds, wfn_fft_grid[0], wfn_fft_grid[1], wfn_fft_grid[2]), dtype='c16')
wfns_folded = np.zeros_like(wfn_ffts) 

for g_idx in range(total_gvecs):
    folded_g = gvecs_folded[g_idx, :]
    wfns_folded[:, folded_g[0], folded_g[1], folded_g[2]] = wfns[:, g_idx]

wfn_ffts = np.fft.ifftn(wfns_folded, axes=(1, 2, 3))
wfn_ffts *= np.sqrt(np.prod(wfn_fft_grid)**2/celvol)   # Times the fft norm factor. 

with open(f'co_band_{bnd_idx}.xsf', 'w') as f:
        xsf.write_xsf(f, [struc], np.abs(wfn_ffts[bnd_idx-1, ...])**2)


# Check normalization. 
wfn_norm = np.sum(np.abs(wfn_ffts[bnd_idx, ...])**2) * celvol/np.prod(wfn_fft_grid)
print(f'Norm of given band index is: {wfn_norm}')
# endregion



# region: Get elph matrix elements.
# First downsample. 
dvscf_downsampled = dvscf[:, ::2, ::2, ::2]    # Downsampled by 2 in real space.
print(dvscf_downsampled.shape)
print(wfn_ffts.shape)

# Then do the umklapp. Normally yeah, but skip for now. 
wfn_umklapped = wfns_folded           # Do the umkapp after this. 

# Finally. 
# Since only one k, q, point. Just n x n' for now. 
nbnds_elph = 20   # Dont need all 500 bands.
elphmat = np.zeros((nat*3, nbnds_elph, nbnds_elph), dtype='c16')
 


# Main workhorse. 

for bnd in range(nbnds_elph):
     for bndp in range(nbnds_elph):
          for mode in range(nat*3):
               
               print(f'Current iter info (n, np, mode): {bnd}, {bndp}, {mode}')

               # region: Calculate electron-phonon matrix entry. 

               # Mat-vec using fft.                     dvscf_uq | \psi_nk >.
               dvscf_wfn = dvscf_downsampled[mode, ...]*wfn_ffts[bndp, ...]                        # Do the real space multiplication. 
               dvscf_wfn_Gspace = np.fft.fftn(dvscf_wfn) * np.sqrt(celvol)/np.prod(wfn_fft_grid)   # Do fft and multiply by a factor. 


               # vec * vec dot product.    < \psi_nk+q | dvscf_uq | \psi_nk >              
               elphmat[mode, bnd, bndp] = np.sum(np.conjugate(wfn_umklapped[bnd, ...]) * dvscf_wfn_Gspace)
               
               # endregion
               
               
# endregion 
x = 1

# region: Test. Print out the format. 
elphmat_double = np.zeros((nat*3, nbnds_elph, nbnds_elph, 2))
elphmat_double[..., 0] = np.real(elphmat)
elphmat_double[..., 1] = np.imag(elphmat)   # order now: (mode, i, j, 2)

elphmat_double = np.einsum('mijd->ijmd', elphmat_double)            # Reorder for format in epw first. 
elphmat_double = elphmat_double.flatten(order='F')

with open('elph.dat', 'w') as f:
     f.write(f'Writing elph matrix.\n')
     f.write(f'Size of the matrix (nbnd, nbnd, nkpts, nmodes, nkpts): {nbnds_elph} {nbnds_elph} 1 {nat*3} 1 \n')
     f.write('Something\n')
     f.write('Something\n')
     f.write('----------------------------------------\n')
     for entry in elphmat_double:
          f.write(str(entry) + '\n')

print('Done writing elph matrix')
 
# endregion. 

# Parallel hdf5 write test. 

# Parallel compilation. hbse test. 