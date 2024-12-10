import matplotlib.pyplot as plt
import math
import csv
import ast
import subprocess
import shutil
import argparse

import os
print(os.environ['PATH'])
os.environ['PATH'] += os.pathsep + r'C:\Program Files\Graphviz\bin'
print(os.environ['PATH'])

import graphviz

from PIL import Image, ImageDraw, ImageFont


# Initialize the parser
parser = argparse.ArgumentParser(description="Process some integers.")

# Adding arguments
parser.add_argument('rep', type=int, help='An integer specifying the repetition number')

# Parse the arguments
args = parser.parse_args()

# Use the arguments
rep = args.rep


print("setting up rep",rep)
# clean up files
folder = 'vis_output/'
# Iterate over each file in the directory
for filename in os.listdir(folder):
    file_path = os.path.join(folder, filename)
    try:
        if os.path.isfile(file_path) or os.path.islink(file_path):
            os.unlink(file_path)  # This removes each file
        elif os.path.isdir(file_path):
            shutil.rmtree(file_path)  # If you want to remove directories, use this
    except Exception as e:
        print(f'Failed to delete {file_path}. Reason: {e}')

# Create cellPredictions directory
os.makedirs('vis_output/cellPredictions', exist_ok=True)  # The exist_ok=True parameter prevents an error if the directory already exists

print("done")
