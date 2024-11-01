from Construct_Runbase import Construct_Runbase
import os

if __name__ == "__main__":
    BASE_DIR = '../ROOT_DIR/'
    ROOT_DIR = '../runbase/'
    NPROC    = 4
    SOURCE_FILE = './SOURCES_SELECTED'

    Construct_Runbase(BASE_DIR,ROOT_DIR,NPROC,SOURCE_FILE)

    os.chdir(ROOT_DIR)
    #os.system('sbatch submit')