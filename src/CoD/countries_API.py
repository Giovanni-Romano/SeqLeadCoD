# WHO Country Data JSON to CSV Converter
# This script converts the JSON from the WHO GHO API to CSV format

import requests
import csv
import json
import pandas as pd

def convert_who_country_data_to_csv():
    """
    Fetches country data from WHO GHO API and converts it to CSV
    """
    # Step 1: Fetch the data
    url = 'https://ghoapi.azureedge.net/api/DIMENSION/COUNTRY/DimensionValues'
    
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad responses
        
        # Parse JSON data
        data = response.json()
        
        # Check if the data has the expected structure
        if 'value' not in data or not isinstance(data['value'], list):
            print('Unexpected data structure')
            return
        
        # Save to CSV
        df = pd.DataFrame(data['value'])
        df.to_csv('data/raw/countries.csv', index=False, sep=';')
        print('CSV file has been saved as data/raw/countries.csv')
    
    except requests.exceptions.RequestException as e:
        print(f'Error fetching data: {e}')
    except Exception as e:
        print(f'Error processing data: {e}')

def main():
    convert_who_country_data_to_csv()

main()