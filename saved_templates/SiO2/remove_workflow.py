import os 

inodes = os.listdir('./')

for inode in inodes:

    # Remove only directory names that start with a number.
    if os.path.isdir(inode) and inode[0].isnumeric():  
        os.system(f'rm -r {inode}')