import os
import subprocess
import platform

print("Current platform: ", current_platform) #print the computer type to terminal

# Check the platform
current_platform = platform.system()

# Run tracking
if current_platform == "Windows":
    subprocess.call('for %i in (*.avi) do sleap-track "%i" --tracking.tracker simple --tracking.similarity centroid -m "centered_instance_model" -m "centroid_model"', shell=True)
elif current_platform == "Darwin":  # macOS
    subprocess.call('for i in *.avi; do sleap-track "$i" --tracking.tracker simple --tracking.similarity centroid -m "centered_instance_model" -m "centroid_model"; done', shell=True)
else:  # Assume Unix-like system
    subprocess.call('for i in *.avi; do sleap-track "$i" --tracking.tracker simple --tracking.similarity centroid -m "centered_instance_model" -m "centroid_model"; done', shell=True)

# Convert tracking to analysis file
if current_platform == "Windows":
    subprocess.call('for %i in (*.slp) do sleap-convert -o "%i.h5" --format analysis "%i"', shell=True)
elif current_platform == "Darwin":  # macOS
    subprocess.call('for i in *.slp; do sleap-convert -o "$i.h5" --format analysis "$i"; done', shell=True)
else:  # Assume Unix-like system
    subprocess.call('for i in *.slp; do sleap-convert -o "$i.h5" --format analysis "$i"; done', shell=True)
