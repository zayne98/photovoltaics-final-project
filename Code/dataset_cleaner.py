# In an effort to get semi-accurate data for irradiance across the US, 
#   I downloaded datasets from https://www.renewables.ninja/.
# This site holds a different dataset for each state in the US, containing
#   weather info including precipitation, temperature, irradiance, snowfall
#   snow mass, cloud cover, and air_density.  As far as I can tell, these 
#   datapoints are averaged across the entire state.  The datasets contain
#   these averages for every day of the year from 1980 through 2019.
# As they are, these datasets are far too massive to be helpful.  Instead,
#   I used this code to process each dataset and strip a year's worth of data
#   instead of 40 years.

# I wish I could have used an online API to regularly ping for just exactly the
#   info I needed.  Sadly, http://www.soda-pro.com/ provides just this feature
#   in an easy to use API, but the feature requires a paid subscription.


# Goal:
#   Only want Temperature and Irradiance_surface, and Air_density
#   Choose 2019 only
#   Combine with other states into single dataset

from os import stat
import pandas as pd
import numpy as np
import glob, os

# Pick one year to use data from
STANDARD_YEAR = 2018

# Path to folder containing the state datasets downloaded from renewables.ninja
FOLDER_PATH = "../state-data/"


# Clean a single state's sheet
def clean_sheet(sheet_name):
    # Read in the sheet
    df = pd.read_csv(sheet_name, header = 2, dtype={1: 'str', 2: float, 
        3: float, 4: float, 5: float,
        6: float, 7: float, 8: float, 
        9: float})

    # We only care about temperature and surface irradiance (and maybe air 
    #   density?). Drop other columns to speed up computations
    df.drop(['precipitation', 'irradiance_toa', 'snowfall', 'snow_mass', 'cloud_cover'], axis=1, inplace=True)
    # print(df.shape)

    # Convert time-stamps to datetime types
    df['time'] =  pd.to_datetime(df['time'])

    # Choose a specific year's worth of data
    df = df[df['time'].dt.year == STANDARD_YEAR]

    df = df.reset_index(drop=True)

    # Tag with state name in a new column
    state_name = sheet_name.split('ninja_weather_country_US.')[1].split("_merra-2")[0]
    print(state_name)
    df["state"] = state_name  
    return df  



# Find all state datasets in the given folder
# I downloaded datasets for all 50 states plus DC
sheet_names = []
os.chdir(FOLDER_PATH)
for file in glob.glob("*.csv"):
    sheet_names.append(file)

combined_df = pd.DataFrame()

for sheet in sheet_names:
    state_df = clean_sheet(sheet)
    combined_df = combined_df.append(state_df)


# print(combined_df.head())
print(combined_df.shape)
# Expect 8760 * 51 rows, and 5 columns

# Remove existing csv and download as new one
os.remove('../Code/irradiance.csv')
combined_df.to_csv('../Code/irradiance.csv', index=False)








