import os,sys
import shutil
sys.path.append('../')
from Utils.Utils import replace_line_in_file

def construct_runbase(BASE_DIR,ROOT_DIR,NPROC,SOURCE_FILE):
    os.system('rm -r %s*' % ROOT_DIR)
    
    shutil.copytree(f'{BASE_DIR}/bin/'         , f'{ROOT_DIR}/bin')
    shutil.copytree(f'{BASE_DIR}/DATA'         , f'{ROOT_DIR}/DATA')
    shutil.copytree(f'{BASE_DIR}/DATABASES_MPI', f'{ROOT_DIR}/DATABASES_MPI')
    shutil.copytree(f'{BASE_DIR}/OUTPUT_FILES' , f'{ROOT_DIR}/OUPUT_FILES')
    
    with open(SOURCE_FILE, 'r') as f:
        info = [line.strip().split()[:3] for line in f]
    nruns = len(info)
    print(len(info))
    
    ### replace the parameters in Par_file
    ### broadcast the mesh and model 
    replace_line_in_file(f'{ROOT_DIR}/DATA/Par_file','NUMBER_OF_SIMULTANEOUS_RUNS     =',f'NUMBER_OF_SIMULTANEOUS_RUNS     = {nruns}')
    replace_line_in_file(f'{ROOT_DIR}/DATA/Par_file','BROADCAST_SAME_MESH_AND_MODEL   =','BROADCAST_SAME_MESH_AND_MODEL   = .true.')
    ### change the simulation type and save the forward wave field 
    replace_line_in_file(f'{ROOT_DIR}/DATA/Par_file','SIMULATION_TYPE                 =','SIMULATION_TYPE                 = 1')
    replace_line_in_file(f'{ROOT_DIR}/DATA/Par_file','SAVE_FORWARD                    =','SAVE_FORWARD                    = .false.')
    
    ### turn the att flag on in our example
    #replace_line_in_file(f'{ROOT_DIR}/DATA/Par_file','ATTENUATION                     =','ATTENUATION                     = .true.')
    
    ### turn the following flags on to save ascii format seis (only format valid in adjoint runs)
    replace_line_in_file(f'{ROOT_DIR}/DATA/Par_file','USE_BINARY_FOR_SEISMOGRAMS      = ','USE_BINARY_FOR_SEISMOGRAMS      = .false.')
    replace_line_in_file(f'{ROOT_DIR}/DATA/Par_file','SU_FORMAT                       = ','SU_FORMAT                       = .false.')
    replace_line_in_file(f'{ROOT_DIR}/DATA/Par_file','ASDF_FORMAT                     = ','ASDF_FORMAT                     = .false.')

    ### displacement is the most usual format in cases (change it if needed)
    replace_line_in_file(f'{ROOT_DIR}/DATA/Par_file','SAVE_SEISMOGRAMS_DISPLACEMENT   =','SAVE_SEISMOGRAMS_DISPLACEMENT   = .true.')
    replace_line_in_file(f'{ROOT_DIR}/DATA/Par_file','SAVE_SEISMOGRAMS_VELOCITY       =','SAVE_SEISMOGRAMS_VELOCITY       = .false.')
    replace_line_in_file(f'{ROOT_DIR}/DATA/Par_file','SAVE_SEISMOGRAMS_ACCELERATION   =','SAVE_SEISMOGRAMS_ACCELERATION   = .false.')
    replace_line_in_file(f'{ROOT_DIR}/DATA/Par_file','SAVE_SEISMOGRAMS_PRESSURE       =','SAVE_SEISMOGRAMS_PRESSURE       = .false.')



    for i, (SOURCE, x, y) in enumerate(info):
        i += 1
        id = str(i).zfill(4)
        source = SOURCE.split('S')[1].zfill(3)
        ### make the simulatuon dirs
        os.mkdir(f'{ROOT_DIR}/run{id}')
        os.mkdir(f'{ROOT_DIR}/run{id}/DATA')
        os.mkdir(f'{ROOT_DIR}/run{id}/DATABASES_MPI')
        os.mkdir(f'{ROOT_DIR}/run{id}/OUTPUT_FILES')
    
        if i == 1:
            os.system('cp %s %s' % (f'{BASE_DIR}/DATABASES_MPI/*', f'{ROOT_DIR}/run{id}/DATABASES_MPI'))
        os.system('cp %s %s' % (f'{BASE_DIR}/OUTPUT_FILES/*', f'{ROOT_DIR}/run{id}/OUTPUT_FILES'))

        ### copy the files needed to dir DATA
        shutil.copy(f'{ROOT_DIR}/DATA/STATIONS' , f'{ROOT_DIR}/run{id}/DATA')
        shutil.copy(f'{ROOT_DIR}/DATA/FORCESOLUTION' , f'{ROOT_DIR}/run{id}/DATA')
        shutil.copy(f'{ROOT_DIR}/DATA/Par_file' , f'{ROOT_DIR}/run{id}/DATA')
        replace_line_in_file(f'{ROOT_DIR}/run{id}/DATA/FORCESOLUTION','FORCE  ',f'FORCE  {source}')
        replace_line_in_file(f'{ROOT_DIR}/run{id}/DATA/FORCESOLUTION','latorUTM:       ',f'latorUTM:       {x}')
        replace_line_in_file(f'{ROOT_DIR}/run{id}/DATA/FORCESOLUTION','longorUTM:      ',f'longorUTM:      {y}')

    os.system('rm %s'% (f'{ROOT_DIR}/DATA/Par_file'))
    

    ### generate the submit file 
    N = NPROC * i
    shutil.copy(f'{BASE_DIR}/submit' , f'{ROOT_DIR}/')
    
    replace_line_in_file(f'{ROOT_DIR}/submit','./bin/xspecfem3D',f'mpirun -np {N} ./bin/xspecfem3D')
    replace_line_in_file(f'{ROOT_DIR}/submit','#SBATCH -n',f'#SBATCH -n {N}')
    replace_line_in_file(f'{ROOT_DIR}/submit','#SBATCH -N',f'#SBATCH -N {i}')

if __name__ == "__main__":
    ### This code is disigned for Specfem 3D Cartesian Version 4.1.0, use it after you have completed
    ### the setting up of your model.  This code is used to set up the parallel computing with the flag 
    ### BROADCAST_SAME_MESH_AND_MODEL on.  
    ### Added by Chao, Li  March 11,2024

    BASE_DIR = '../../XSH_solver/'
    ROOT_DIR = '../../runbase/'
    NPROC    = 4
    SOURCE_FILE = '/public/home/acapp7905d/li_chao/specfem3d-kunshan/Test_li_chao/Tools4SEM3D_MultiSever/files/SOURCES_SELECTED_INTERATION1'

    construct_runbase(BASE_DIR,ROOT_DIR,NPROC,SOURCE_FILE)
