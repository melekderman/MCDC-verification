import h5py

file_path = "output.h5"

with h5py.File(file_path, 'r') as h5_file:
    def hdf5_group(group, indent=0):
        for key in group:
            item = group[key]
            if isinstance(item, h5py.Group):
                print(" " * indent + f"Group: {key}")
                hdf5_group(item, indent + 4)
            elif isinstance(item, h5py.Dataset):
                print(" " * indent + f"Dataset: {key}, Shape: {item.shape}")
                if item.shape == ():  
                    print(" " * (indent + 4) + f"Scalar value: {item[()]}")
                else:  
                    print(" " * (indent + 4) + f"Data: {item[:]}")
            else:
                print(" " * indent + f"Unknown: {key}")

    hdf5_group(h5_file)
