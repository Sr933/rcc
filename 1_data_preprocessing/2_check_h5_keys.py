import h5py
import os
main_folder="/data/sr933/scRCC validation/GSE159115_RAW"

for folder in os.listdir(main_folder):
    if ".h5" in folder:
        folder_path=os.path.join(main_folder, folder)

        with h5py.File(folder_path, 'r') as hdf:
            print(list(hdf.keys()))
