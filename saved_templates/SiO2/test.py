from ase import atoms, io


sio2 = io.read('SiO2.cif')

print(type(sio2))

print(type(sio2.cell.cellpar().tolist()))

print(sio2.get_pbc())

print(sio2.get_positions().tolist())