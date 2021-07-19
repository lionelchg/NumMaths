import os

def create_dir(folder):
    if not os.path.exists(folder): 
        os.makedirs(folder)