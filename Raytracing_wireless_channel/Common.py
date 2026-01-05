from sionna.rt import Transmitter,PathSolver
import numpy as np
import scipy.io
import os

def CalculateMeanPathPower(scene,APNum,UEPositions):
    
    old_txs = list(scene.transmitters.keys())
    for tx_name in old_txs:
        scene.remove(tx_name)
    
    UENum = UEPositions.shape[1]
    for UEidx in range(UENum):
        scene.add(Transmitter(name=f"UE{UEidx}", position=UEPositions[:,UEidx],color=(1.0, 0.0, 0.0), display_radius=0.2))
     
    LoSPathPowers = []
    SingleReflectionPathPowers = []
    DoubleReflectionPathPowers = []
    SingleScatteringPathPowers = []
    SingleDiffractionPathPowers = []
    SingleRefractionPathPowers = []
    ScatteringReflectionPathPowers = []
    ReflectionScatteringPathPowers = []
    DoubleScatteringPathPowers = []
    
    p_solver = PathSolver()                    
    paths = p_solver(scene=scene, 
                    max_depth=2,
                    max_num_paths_per_src=int(5e5),
                    samples_per_src=int(5e5),
                    synthetic_array=True, 
                    los=True, 
                    specular_reflection=True,
                    diffuse_reflection=True, 
                    refraction=True, 
                    diffraction=True,
                    seed=int(42))
    coefficients_real = paths.a[0].numpy()
    coefficients_imaginary = paths.a[1].numpy()
    coefficients = coefficients_real + 1j*coefficients_imaginary
    PathInteractionTypes = paths.interactions.numpy()
    PathNum = PathInteractionTypes.shape[3]
    for UEidx in range(UENum):
        for APidx in range(APNum):
            for Pathidx in range(PathNum):
                PathPower = np.abs(coefficients[APidx,0,UEidx,0,Pathidx])**2
                PathType = PathInteractionTypes[:,APidx,UEidx,Pathidx]
                # DIFFRACTION = 8, REFRACTION = 4, DIFFUSE = 2, SPECULAR = 1, NONE = 0
                if (PathType==[0,0]).all() and PathPower>0:
                    LoSPathPowers.append(PathPower)
                elif (PathType==[1,0]).all():
                    SingleReflectionPathPowers.append(PathPower)
                elif (PathType==[1,1]).all():
                    DoubleReflectionPathPowers.append(PathPower)    
                elif (PathType==[2,0]).all():
                    SingleScatteringPathPowers.append(PathPower)
                elif (PathType==[8,0]).all():
                    SingleDiffractionPathPowers.append(PathPower)
                elif (PathType==[4,0]).all():
                    SingleRefractionPathPowers.append(PathPower)
                elif (PathType==[1,2]).all():
                    ReflectionScatteringPathPowers.append(PathPower)
                elif (PathType==[2,1]).all():
                    ScatteringReflectionPathPowers.append(PathPower)
                elif (PathType==[2,2]).all():
                    DoubleScatteringPathPowers.append(PathPower)                            
                       
    return [LoSPathPowers,SingleReflectionPathPowers,DoubleReflectionPathPowers,SingleScatteringPathPowers,SingleDiffractionPathPowers,SingleRefractionPathPowers,ScatteringReflectionPathPowers,ReflectionScatteringPathPowers,DoubleScatteringPathPowers]

def GenCIRSamples(scene, UEPositions, PickedPathNum, APPositions, SpaceBound, filename): 

    LightSpeed = 3e8
    APNum = APPositions.shape[1]
    
    old_txs = list(scene.transmitters.keys())
    for tx_name in old_txs:
        scene.remove(tx_name)
        
    UENum = UEPositions.shape[1]
    
    for UEidx in range(UENum):
        pos = UEPositions[:, UEidx]
        scene.add(Transmitter(name=f"UE{UEidx}", position=pos))
            
    p_solver = PathSolver()                    
    paths = p_solver(scene=scene, 
                    max_depth=2, 
                    max_num_paths_per_src=int(5e5),
                    samples_per_src=int(5e5),       
                    synthetic_array=True, 
                    los=True, 
                    specular_reflection=True,
                    diffuse_reflection=True, 
                    refraction=True, 
                    diffraction=True,
                    seed=int(42))
    
    taus = paths.tau.numpy()
    coeffs_real = paths.a[0].numpy()
    coeffs_imag = paths.a[1].numpy()
    coeffs_complex = coeffs_real + 1j * coeffs_imag
    inter_types = paths.interactions.numpy() 
    inter_objects = paths.objects.numpy()
    inter_points = paths.vertices.numpy()
    
    OutputData = []
    Room_struct = {
                    'SpaceBound': SpaceBound,
                    'APPositions': APPositions
                } 
    OutputData.append(Room_struct)
    
    for UEidx in range(UENum):   
        cir_cell_array = np.empty((APNum, 1), dtype=object)
        
        for APidx in range(APNum):            
            
            raw_delays = taus[APidx, UEidx, :]*LightSpeed
            raw_coeffs = coeffs_complex[APidx, 0, UEidx, 0, :]
            raw_types  = inter_types[:, APidx, UEidx, :] 
            raw_objects = inter_objects[:, APidx, UEidx, :]
            raw_points = inter_points[:, APidx, UEidx, :, :]
            
            raw_powers = np.abs(raw_coeffs)**2
            
            nonzero_indices = np.flatnonzero(raw_powers > 0)
            sorted_relative_indices = np.argsort(-raw_powers[nonzero_indices])
            valid_indices = nonzero_indices[sorted_relative_indices[0:PickedPathNum]]#要注意先取出功率非零的有效径，再取PickedPathNum条功率最大的出来            
            
            ap_paths = np.zeros((0, 13))
            if len(valid_indices) > 0:
                ap_paths = np.zeros((len(valid_indices),13))
                
                curr_delays = raw_delays[valid_indices]
                curr_coeffs = raw_coeffs[valid_indices]
                curr_types  = raw_types[:,valid_indices]
                curr_objects = raw_objects[:,valid_indices]   
                curr_points = raw_points[:,valid_indices,:]
                
                for k in range(len(valid_indices)):
                    ap_paths[k,:] = np.r_[curr_delays[k], np.real(curr_coeffs[k]), np.imag(curr_coeffs[k]), curr_types[:,k], 
                                     curr_objects[:,k], curr_points[0,k,:], curr_points[1,k,:]]

            cir_cell_array[APidx, 0] = ap_paths
                
        ue_struct = {
            'UE_Pos': UEPositions[:, UEidx],
            'CIR_by_AP': cir_cell_array
        }                    
        OutputData.append(ue_struct)
        
    base_dir = os.path.dirname(os.path.abspath(__file__))
    save_dir = os.path.join(base_dir, 'FigureDatas')
    os.makedirs(save_dir, exist_ok=True) 
    save_path = os.path.join(save_dir, f"{filename}.mat")
    scipy.io.savemat(save_path, {filename: OutputData})
    print(f"Data successfully saved to {filename}.mat")