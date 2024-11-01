from obspy import Trace, UTCDateTime
import numpy as np

def replace_line_in_file(file_in,search_string, new_content):
    with open(file_in, "r") as file:
        lines = file.readlines()
        
    with open(file_in, "w") as file:
        for line in lines:
            if search_string in line:
                line = new_content + '\n'
            file.write(line)

### calculate the signal/noise ratio
def __rms(arr):
    return np.sqrt(np.mean(np.square(arr)))

def __SNR__(t, data, t1, t2, noi_len):
    sig_win = np.where((t > t1) & (t < t2) )
    noi_win = np.where((t > t2+1) & (t < t2+1+noi_len))
    SNR =  __rms(data[sig_win])/__rms(data[noi_win])
    return SNR

### Tools to construct the runbase
def replace_line_in_file(file_in,search_string, new_content):
    with open(file_in, "r") as file:
        lines = file.readlines()
        
    with open(file_in, "w") as file:
        for line in lines:
            if search_string in line:
                line = new_content + '\n'
            file.write(line)

def find_line_in_file(file_in,search_string):
    with open(file_in, "r") as file:
        lines = file.readlines()
        
    with open(file_in, "r") as file:
        for line in lines:
            if search_string in line:
                dt = line.strip().split('=')[1]
    return dt

### I/O of the data 
def write_dat(file, t, data):
    f1 = open(file, 'w')
    for i, T in enumerate(t):
        f1.write(str(T)+' '+str(data[i])+'\n')
    f1.close
    #data = np.array(data,dtype=np.float32)
    #data.tofile(file)

def read_dat(file):
    data = np.loadtxt(file)
    dt   = data[3,0] - data[2,0]
    data = data[3:-1, 1] + data[3:-1, 2]
    egf  = np.gradient(data, dt)

    return egf

def read_semd(file):
    data = np.loadtxt(file)
    t    = data[:,0]
    syn  = data[:,1]
    return t,syn

### Read the window par for cal adj src
def get_info(config):
    with open(config.window_file, 'r') as f1:
        periods = [[float(x) for x in item.split()[:2]] for item in f1.readline().split('\n')[0].split(';')[:3]]
        vel_win = [[float(x) for x in item.split()[:2]] for item in f1.readline().split('\n')[0].split(';')[:3]]
        tshift  = [[float(x) for x in item.split()[:2]] for item in f1.readline().split('\n')[0].split(';')[:3]]

    info = {'periods': periods, 'vel_win': vel_win, 'tshift': tshift}
    return info

### Read the parameters from par_file
def create_dict(filtered_lines):
    result_dict = {}

    for line in filtered_lines:
        columns = line.split()
        if len(columns) >= 3:  
            key = columns[0]   
            value = columns[2]  
            result_dict[key] = value  

    return result_dict

def read_par(file):
    with open(file) as f1:
        lines = f1.readlines()

    filtered_lines = [line.strip() for line in lines if not line.strip().startswith('#') and line.strip()]

    par = create_dict(filtered_lines)

    class item:
        def __init__(self):
            self.sources          = par["sources"]
            self.obs_path         = par["obs_path"]
            self.syn_path         = par["syn_path"]
            self.window_file      = par["window_file"]

            self.NPROC            = par["NPROC"]

            self.adjsrc_type      = par["adjsrc_type"]

            self.base_path        = par["base_path"]

            self.half_duration    = np.float(par["half_duration"]   )
            self.noi_len          = np.float(par["noi_len"]         )

            self.taper_percentage = np.float(par["taper_percentage"])
            self.taper_type       = par["taper_type"]      
            self.measure_type     = par["measure_type"]    
            
            self.dt_sigma_min     = np.float(par["dt_sigma_min"]    )
            self.dlna_sigma_min   = np.float(par["dlna_sigma_min"]  )

            self.snr              = np.float(par["snr"]             )
            self.tshift           = np.float(par["tshift"]          )
    config = item()

    return config