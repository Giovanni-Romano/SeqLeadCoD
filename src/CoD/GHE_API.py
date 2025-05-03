import requests
import itertools
import pandas as pd
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import urllib.parse

# Function to fetch country codes
def get_country_codes():
    url = "https://xmart-api-public.who.int/DEX_CMS/GHE_FULL?$apply=groupby((DIM_COUNTRY_CODE))"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        return [item["DIM_COUNTRY_CODE"] for item in data["value"]]
    except Exception as e:
        print(f"Error fetching country codes: {str(e)}")
        return []

# Function to fetch age group codes
def get_age_group_codes():
    url = "https://xmart-api-public.who.int/DEX_CMS/GHE_FULL?$apply=groupby((DIM_AGEGROUP_CODE))"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        return [item["DIM_AGEGROUP_CODE"] for item in data["value"]]
    except Exception as e:
        print(f"Error fetching age group codes: {str(e)}")
        return []

# Fetch country codes and age group codes
country_codes = get_country_codes()
age_groups = get_age_group_codes()

# Print the fetched codes for verification
print(f"Fetched {len(country_codes)} country codes: {country_codes}")
print(f"Fetched {len(age_groups)} age group codes: {age_groups}")

sexes = ['MALE', 'FEMALE']
years = list(range(2000, 2022))
base_url = "https://xmart-api-public.who.int/DEX_CMS/GHE_FULL"

# 2. Prepare all combinations
combos = list(itertools.product(country_codes, years, sexes, age_groups))

# 3. Function to fetch data
def fetch_data(combo):
    country, year, sex, age = combo
    
    try:
        # Create filter parameter
        filter_param = f"FLAG_RANKABLE eq 1 and DIM_COUNTRY_CODE eq '{country}' and DIM_SEX_CODE eq '{sex}' and DIM_AGEGROUP_CODE eq '{age}' and DIM_YEAR_CODE eq {year}"
        
        # Fixed parameter order and selection
        params = {
        "$filter": filter_param,
        "$orderBy": "VAL_DTHS_RATE100K_NUMERIC desc",
        "$select": "DIM_COUNTRY_CODE, DIM_SEX_CODE, DIM_AGEGROUP_CODE, DIM_YEAR_CODE, DIM_GHECAUSE_CODE,DIM_GHECAUSE_TITLE, VAL_DTHS_RATE100K_NUMERIC",
        "$top": 1
        }
        
        encoded_params = urllib.parse.urlencode(params)
        final_url = f"{base_url}?{encoded_params}"
        
        response = requests.get(final_url)
        response.raise_for_status()  # This will trigger the exception for bad responses
        data = response.json()
        
        # Brief pause to avoid server throttling
        time.sleep(0.2)
        
        if 'value' in data and data['value']:
            return data['value']
        else:
            return []
    
    except Exception as e:
        # Return empty list and log the error separately
        print(f"Error with {country}, {year}, {sex}, {age}: {str(e)}")
        raise


# 4. Set up result tracking
results = []
failed_combos = []
exception_combos = []
successful_count = 0
failed_count = 0
exception_count = 0


# 5. Use semaphore to limit concurrent connections
max_workers = 10  # Reduced from 10 to prevent overwhelming the API

with ThreadPoolExecutor(max_workers=max_workers) as executor:
    # Submit all jobs
    futures = {executor.submit(fetch_data, combo): combo for combo in combos}
    
    # Process results as they complete
    for future in tqdm(as_completed(futures), total=len(futures), mininterval=5, desc="Fetching WHO GHE data"):
        combo = futures[future]
        try:
            result = future.result()
            if result:  # If we got data
                results.extend(result)
                successful_count += 1
            else:
                failed_count += 1
                failed_combos.append(combo)
                print(f"Future failed for {combo}: {str(e)}")
        except Exception as e:
            exception_combos.append(combo)
            print(f"Future failed for {combo}: {str(e)}")
            exception_count += 1

# 6. Save data
if results:
    df = pd.DataFrame(results)
    df.to_csv("data/WHO_GHE_Top1_Parallel.csv", index=False)
    print(f"Saved {len(df)} rows to WHO_GHE_Top1_Parallel.csv.")
else:
    print("No results collected.")

print(f"Statistics: {successful_count} successful queries, {failed_count} failed queries")