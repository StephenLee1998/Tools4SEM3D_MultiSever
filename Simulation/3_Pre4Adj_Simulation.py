import os
import shutil
from Utils.Utils import read_par, replace_line_in_file

def pre4adj_simulation(config):
    
    with open(config.sources, 'r') as f:
        info = [line.strip().split()[:3] for line in f]

    for i, (SOURCE, x, y) in enumerate(info):
        i += 1
        id = str(i).zfill(4)
        source = SOURCE.split('S')[1].zfill(3)
        
        print(source,id)
        adj_sta = './scratch/traces/S'+source+'/STATIONS_ADJOINT'
        work_sta = config.syn_path+'/run'+id+'/DATA/STATIONS_ADJOINT'
        adj_dir = './scratch/traces/S'+source+'/adj'
        work_dir = config.syn_path+'/run'+id+'/SEM'

        par_file = config.syn_path+'/run'+id+'/DATA/Par_file'

        replace_line_in_file(par_file,'SIMULATION_TYPE                 =','SIMULATION_TYPE                 = 3')
        replace_line_in_file(par_file,'SAVE_FORWARD                    =','SAVE_FORWARD                    = .false.')
        #replace_line_in_file(par_file,'ATTENUATION                     =','ATTENUATION                     = .false.')

        shutil.copy2(adj_sta,work_sta)

        if os.path.exists(work_dir):
            shutil.rmtree(work_dir)
            shutil.copytree(adj_dir,work_dir)
        else:
            shutil.copytree(adj_dir,work_dir)



if __name__ == "__main__":
    ### This code is disigned for Specfem 3D Cartesian Version 4.1.0, use it after you have calculated
    ### the adjoint sources. 

    config = read_par('./par_file')

    pre4adj_simulation(config)
