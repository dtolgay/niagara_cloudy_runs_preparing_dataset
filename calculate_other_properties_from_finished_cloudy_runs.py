import sys
sys.path.append("/scratch/m/murray/dtolgay")


import numpy as np 
import pandas as pd 

from tools import constants

def main(base_file_dir):
    centers = read_centers_file(base_file_dir=base_file_dir)

    properties = []
    for row, center in centers.iterrows():

        metallicity_center = 10**center['log_metallicity'] 

        try: 
            densities = read_ovr_file(
                base_file_dir=base_file_dir,
                center=center
                )
            
            average_fh2 = calculate_fh2(densities)
            average_fCO = calculate_fCO(
                densities=densities, 
                metallicity=metallicity_center
                )

            properties.append((average_fh2, average_fCO))
        except Exception as e: 
            # If exception occurs return NaN
            properties.append((np.nan, np.nan))
            print(f"Exception occured for run {center}: \n{e}")

        if row % 1e4 == 1: 
            print(f"{row} finished. Left {len(centers) - row}")


    centers[['fh2', 'fCO']] = np.array(properties)

    print(centers)

    write_to_a_file(
        df = centers,
        base_file_dir = base_file_dir,
        file_name = "other_properties.csv"
    )

    return 0

def read_centers_file(base_file_dir):

    # Get the file path
    centers_file_path = f"{base_file_dir}/centers.txt"
    centers_train = np.loadtxt(fname=centers_file_path) 

    print(f"Using {centers_file_path} as a center.txt file")

    centers_train_df = pd.DataFrame(
        centers_train,
        columns=[
            "log_metallicity",
            "log_hden",
            "log_turbulence",
            "log_isrf",
            "log_radius",
        ],
    )

    return centers_train_df

def read_ovr_file(base_file_dir, center):

    fdir = f"hden{center['log_hden']:.5f}_metallicity{center['log_metallicity']:.5f}_turbulence{center['log_turbulence']:.5f}_isrf{center['log_isrf']:.5f}_radius{center['log_radius']:.5f}"

    # Read files 
    densities = pd.read_csv(f"{base_file_dir}/{fdir}/{fdir}.ovr", delim_whitespace=True)
    densities.rename(columns={'#depth': 'depth'}, inplace=True)

    # Calculate the column density 
    densities['column_density'] = \
    densities['depth'] * 10**center['log_hden'] * constants.proton_mass * constants.mu_h * constants.kg2g # gr / cm^2    

    return densities

def calculate_fh2(densities):

    # Calculate total h2 column density  
    sigma_h2 = np.zeros(len(densities)) * np.nan
    sigma_h = np.zeros(len(densities)) * np.nan
    depth_difference = np.zeros(len(densities)) * np.nan

    for row, density in densities.iterrows():
        if row == 0:
            sigma_h2[row] = density['2H_2/H'] * density['hden'] * constants.proton_mass * constants.kg2g * density['depth']
            sigma_h[row] = density['hden'] * constants.proton_mass * constants.kg2g * density['depth']
            depth_difference[row] = density['depth']
        else: 
            sigma_h2[row] = density['2H_2/H'] * density['hden'] * constants.proton_mass * constants.kg2g * (density['depth'] - previous_row['depth'])
            sigma_h[row] = density['hden'] * constants.proton_mass * constants.kg2g * (density['depth'] - previous_row['depth'])
            depth_difference[row] = density['depth'] - previous_row['depth']
        
        # Store the previous row for depth calculation
        previous_row = density
        
    fh2 = sigma_h2 / sigma_h    
    averaged_fh2 = sum(fh2 * sigma_h) / sum(sigma_h)

    return averaged_fh2

def calculate_fCO(densities, metallicity):
    
    C_over_H_mass_ratio = 12 / 1 
    C_over_H_number_ratio = 2.51e-4 # ~ C/H /home/m/murray/dtolgay/cloudy/c23.00/data/abundances/ISM.abn 
    CO_over_C_mass_ratio = 28 / 12
    
    C_mass_ratio = C_over_H_number_ratio * C_over_H_mass_ratio 

    
    # Set hden
    hden = densities.iloc[0]['hden'] * constants.proton_mass * constants.kg2g # hden is constant in the run [gr/cm^3]
    C_density = hden * C_mass_ratio * metallicity # Scaling according to the metallicity [g/cm^3] 
    
    # Initiate the arrays 
    sigma_C = np.zeros(len(densities)) * np.nan
    sigma_h = np.zeros(len(densities)) * np.nan
    sigma_CO = np.zeros(len(densities)) * np.nan
    sigma_h2 = np.zeros(len(densities)) * np.nan
    

    # Calculate CO abundance
    for row, density in densities.iterrows():
        if row == 0:
            sigma_C[row] = C_density * density['depth']
            
            sigma_CO[row] = sigma_C[row] * density['CO/C'] * CO_over_C_mass_ratio
            
            sigma_h[row] = hden * density['depth']
            
            sigma_h2[row] = density['2H_2/H'] * hden * density['depth']            
        
        else: 
            sigma_C[row] = C_density * (density['depth'] - previous_row['depth'])
            
            sigma_CO[row] = sigma_C[row] * density['CO/C'] * CO_over_C_mass_ratio
            
            sigma_h[row] = hden * (density['depth'] - previous_row['depth'])
            
            sigma_h2[row] = density['2H_2/H'] * hden * (density['depth'] - previous_row['depth'])           

        # Store the previous row for depth calculation
        previous_row = density

    f_CO = sigma_CO / sigma_h
    averaged_f_CO = sum(f_CO * sigma_h) / sum(sigma_h) 

    return averaged_f_CO

def write_to_a_file(df, base_file_dir, file_name):

    fname = f"{base_file_dir}/{file_name}"

    df.to_csv(
        fname,
        index=False,
        )

    print(f"File is written to {fname}")

    # Can be read as: df_read = pd.read_csv(fname) 


    return 0


if __name__ == "__main__":

    # base_file_dir = "/scratch/m/murray/dtolgay/cloudy_runs/z_0/cr_1_CO87_CII_H_O3/cr_1_CO87_CII_H_O3_metallicity_minus2_minus3point5"
    base_file_dir = "/scratch/m/murray/dtolgay/cloudy_runs/z_0/cr_1_CO87_CII_H_O3/cr_1_CO87_CII_H_O3_metallicity_above_minus_2"


    main(base_file_dir)