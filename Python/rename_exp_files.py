import os
import glob
from tkinter import Tk, filedialog

# --- Step 1: Select folder
Tk().withdraw()  # Hide the root window
folder_path = filedialog.askdirectory(title='Select folder')

if not folder_path:
    raise Exception("No folder selected.")

# --- Step 2: Get list of *_1.avi files
avi_files = glob.glob(os.path.join(folder_path, '*_1.avi'))
avi_names = [os.path.basename(f) for f in avi_files]

if not avi_names:
    raise Exception("No matching *_1.avi files found.")

# --- Step 3: Ask user to choose a file name string to replace
print("Select a string to replace:")
for i, name in enumerate(avi_names):
    print(f"{i+1}: {name}")

index = int(input("Enter the number corresponding to the file: ")) - 1
startstring = avi_names[index][:-6]  # Remove '_1.avi'

# --- Step 4: Find all matching files
all_files = os.listdir(folder_path)
target_files = [f for f in all_files if startstring in f]

# --- Step 5: Rename
name_replace = 'C2_survival_40_empty'

for old_name in target_files:
    new_name = old_name.replace(startstring, name_replace)
    old_path = os.path.join(folder_path, old_name)
    new_path = os.path.join(folder_path, new_name)
    os.rename(old_path, new_path)

print("Finished renaming files.")
