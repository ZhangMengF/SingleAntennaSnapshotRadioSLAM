import sionna.rt
from sionna.rt import load_scene,PlanarArray,Receiver,RadioMaterial
import matplotlib.pyplot as plt
import numpy as np
import os
from Common import GenCIRSamples
from Common import CalculateMeanPathPower

# Create scene
print(f"Sionna-rt version: {sionna.rt.__version__}") 
base_dir = os.path.dirname(os.path.abspath(__file__))
scene = load_scene(f"{base_dir}/MeetingRoom/MeetingRoom.xml", merge_shapes=False) 

scene.tx_array = PlanarArray(num_rows=1, 
                            num_cols=1, 
                            vertical_spacing=0.5, 
                            horizontal_spacing=0.5, 
                            pattern="iso", 
                            polarization="VH", 
                            polarization_model="tr38901_2") 
scene.rx_array =scene.tx_array    

scene.frequency = 10e9

HumanSkin = RadioMaterial( 
                    name ="HumanSkin", 
                    relative_permittivity = 31.3, 
                    conductivity = 8.01, 
                    thickness = 0.5,  
                    scattering_coefficient = 1) 
scene.get("Persons").radio_material = HumanSkin

scene.get("Tables").radio_material.scattering_coefficient = 1
scene.get("Tables").radio_material.thickness = 0.1
scene.get("Chairs").radio_material.scattering_coefficient = 1
scene.get("Chairs").radio_material.thickness = 0.1
scene.get("Door").radio_material.scattering_coefficient = 1
scene.get("Door").radio_material.thickness = 0.1
for WallName in ['Wall1','Wall2','Wall3','Wall4','Floor','Ceiling']:
    scene.get(WallName).radio_material.scattering_coefficient = 0
    scene.get(WallName).radio_material.thickness = 1
    # print(scene.get(WallName).object_id)
    # print(scene.get(WallName).position.numpy())

for obj in scene.objects.values(): 
    print(f"{obj.name}: {obj.radio_material}") 

# scene.preview()

## Generate positions of APs and UEs
x_vals = np.arange(-4, 5) 
y_vals = np.r_[-4:-1, 2:5] 
X_grid, Y_grid = np.meshgrid(x_vals, y_vals)
part1_x, part1_y = X_grid.flatten(), Y_grid.flatten()
part2_y = np.arange(-1, 2)
part2_x = np.full_like(part2_y, 4, dtype=float)
all_x = np.concatenate([part1_x, part2_x])
all_y = np.concatenate([part1_y, part2_y])
xy_array = np.vstack((all_x, all_y))
UE_heights = np.full(xy_array.shape[1],0.8)
UEPositions = np.vstack((xy_array,UE_heights))
UENum = UEPositions.shape[1]

Length,Width,Height = 10.0,10.0,5.0
SpaceBound = np.array([Length,Width,Height],dtype=np.float64)
APPositions = np.array([[0,0,5],[0,4,5],[2*np.sqrt(3),2,5],[2*np.sqrt(3),-2,5],[0,-4,5],[-2*np.sqrt(3),-2,5],[-2*np.sqrt(3),2,5]],dtype=np.float64).T 
APNum = APPositions.shape[1]
old_rxs = list(scene.receivers.keys())
for rx_name in old_rxs:
    scene.remove(rx_name)
for APidx in range(APNum):
    scene.add(Receiver(name=f"AP{APidx}", position=APPositions[:,APidx], color=(0.0, 0.0, 1.0), display_radius=0.2))
    
## Generate multipaths
scene.get("Tables").radio_material.scattering_coefficient = 1
scene.get("Chairs").radio_material.scattering_coefficient = 1
scene.get("Door").radio_material.scattering_coefficient = 1
for WallName in ['Wall1','Wall2','Wall3','Wall4','Floor','Ceiling']:
    scene.get(WallName).radio_material.scattering_coefficient = 0.4
scene.get("Persons").radio_material.scattering_coefficient = 1

PickedPathNum = 100
GenCIRSamples(scene, UEPositions, PickedPathNum, APPositions, SpaceBound, "RaytracedPaths")

## Calculate the power distribution of different types of paths
scene.get("Tables").radio_material.scattering_coefficient = 1
scene.get("Chairs").radio_material.scattering_coefficient = 1
scene.get("Door").radio_material.scattering_coefficient = 1
for WallName in ['Wall1','Wall2','Wall3','Wall4','Floor','Ceiling']:
    scene.get(WallName).radio_material.scattering_coefficient = 0.4
scene.get("Persons").radio_material.scattering_coefficient = 1

MultiTypePathPowers = CalculateMeanPathPower(scene,APNum,UEPositions)

PathTypeNum = len(MultiTypePathPowers)
TotalPathNum = 0
TitleList = ['LoS','Single Reflection','Double Reflection',\
    'Single Scattering','Single Diffraction','Single Refraction','Scattering-Reflection',\
    'Reflection-Scattering','Double Scattering']

for PathTypeIdx in range(PathTypeNum):
    PathPowers = np.array(MultiTypePathPowers[PathTypeIdx])
    TotalPathNum += len(PathPowers)
    plt.hist(PathPowers, bins=100)
    plt.title(f'{scene.frequency.numpy().item():.1e}Hz '+TitleList[PathTypeIdx]+f' Path Powers Distribution\n mean:{PathPowers.mean():.2e},std:{PathPowers.std():.2e}, PathNum per UE per AP:{len(PathPowers)/(APNum*UENum):.1e}')
    plt.savefig(f'{base_dir}/FigureDatas/{scene.frequency.numpy().item():.1e}Hz_{TitleList[PathTypeIdx]} PathPowers.png')
    plt.close()
print(f'Average path number per UE per AP: {TotalPathNum/(APNum*UENum):.1e}')
print("Calculation of path power distribution completed.")    