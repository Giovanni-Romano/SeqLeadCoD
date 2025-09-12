import requests
import itertools
import pandas as pd
import time
import urllib.parse
from tqdm import tqdm  # Added import for tqdm

# 1. Define country codes
country_codes = [
    'FIN'  # Changed from AFG to ALB as in your desired URL
]

# 2. Define age group codes
age_groups = [
    'TOTAL'  # Changed from Y0T6 to Y20T24 as in your desired URL
]
# Uncomment to use more age groups
# age_groups = ['Y0T6', 'Y7T27D', 'Y28D364D', 'Y1T4', 'Y5T9', 'Y10T14', 'Y15T19', 'Y20T24',
#              'Y25T29', 'Y30T34', 'Y35T39', 'Y40T44', 'Y45T49', 'Y50T54', 'Y55T59', 'Y60T64',
#              'Y65T69', 'Y70T74', 'Y75T79', 'Y80T84', 'Y85T89', 'Y90T94', 'Y95T99', 'YGE100']

# 3. Sex codes (only MALE and FEMALE)
sexes = ['FEMALE']
# Uncomment to include female
# sexes = ['MALE', 'FEMALE']

# 4. Years (choose your range)
years = [2020]  # Changed from 2000 to 2015 as in your desired URL

# 5. Base URL
base_url = "https://xmart-api-public.who.int/DEX_CMS/GHE_FULL"
# base_url = "https://xmart-api-public-uat.who.int/DEX_CMS/GHE_COD"

# 6. Collect all results
all_results = []

combos = list(itertools.product(country_codes, years, sexes, age_groups))

country = combos[0][0]  # SLB
year = combos[0][1]  # 2004
sex = combos[0][2]  # FEMA
age = combos[0][3]  # Y20T24

# 7. Loop over combinations
for country, year, sex, age in tqdm(combos, desc="Fetching WHO GHE data"):
    # Create filter parameter
    filter_param = f"FLAG_RANKABLE eq 1 and DIM_COUNTRY_CODE eq '{country}' and DIM_SEX_CODE eq '{sex}' and DIM_AGEGROUP_CODE eq '{age}' and DIM_YEAR_CODE eq {year}"
    
    # Changed parameter order and selection to match your desired URL
    params = {
        "$filter": filter_param,
        "$orderBy": "VAL_DTHS_RATE100K_NUMERIC desc",
        "$select": "DIM_COUNTRY_CODE, DIM_SEX_CODE, DIM_AGEGROUP_CODE, DIM_YEAR_CODE, DIM_GHECAUSE_CODE,DIM_GHECAUSE_TITLE, VAL_DTHS_RATE100K_NUMERIC",
        "$top": 1
    }
    
    encoded_params = urllib.parse.urlencode(params)
    final_url = f"{base_url}?{encoded_params}"
    
    print(f"Generated URL: {final_url}")  # Print the URL for verification
    
    try:
        response = requests.get(final_url)
        response.raise_for_status()
        data = response.json()
        
        if 'value' in data:
            all_results.extend(data['value'])
        else:
            print(f"No data for {country}, {year}, {sex}, {age}")
    
    except Exception as e:
        print(f"Error with {country}, {year}, {sex}, {age}: {e}")
    
    time.sleep(0.2)  # Sleep to avoid hammering the server

# 8. Save to CSV
if all_results:
    df = pd.DataFrame(all_results)
    df.to_csv("WHO_GHE_Top1.csv", index=False)
    print("Saved all Top 1 data to WHO_GHE_Top1.csv")
else:
    print("No results found. Check the API response.")
# %%
