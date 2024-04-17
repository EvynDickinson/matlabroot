
import os
import shutil
import time

# Set up directories
curr_folder = r'G:\My Drive\Jeanne Lab\DATA'
target_dir = r'D:\Processed Data'

# Get list of data folders
folder_names = [name for name in os.listdir(curr_folder) if os.path.isdir(os.path.join(curr_folder, name)) and len(name) == 10]

# Copy the folders
#for folder_name in folder_names:
#    letter_list = ['A', 'B', 'C', 'D']
#    for letter in letter_list:
#        specific_folder = os.path.join(curr_folder, folder_name, f'Arena {letter}')
#        # Create target directory if it doesn't exist
#        os.makedirs(target_folder, exist_ok=True)
#        shutil.copytree(specific_folder, target_folder)
#    print(f"Copied {folder_name}")


# Select only the first element of folder_names list
folder_name = folder_names[0]

start_time = time.time()

letter_list = ['A', 'B', 'C', 'D']
for letter in letter_list:
    specific_folder = os.path.join(curr_folder, folder_name, f'Arena {letter}')
    target_folder = os.path.join(target_dir, folder_name, f'Arena {letter}')
    # Create target directory if it doesn't exist
    os.makedirs(target_folder, exist_ok=True)
    shutil.copytree(specific_folder, target_folder)

end_time = time.time()

print(f"Copied {folder_name}")
print(f"Time taken: {end_time - start_time} seconds")
