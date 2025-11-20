
#Status codes
#200 OK: Request successful, data returned.
#201 Created: New resource created.
#204 No Content: Success but no data returned.
#400 Bad Request: Invalid request.
#401 Unauthorized: Missing or invalid API key.
#404 Not Found
#500 Internal Server Error: Server encountered an error.

#put = adds or amends
#get = retrieve data (no body)
#post = send data (must have a body - data you're sending)
#del = deletes data 

import requests
import json
import html
import re

def sanitise_json(data):
    """
    Recursively sanitize all strings inside the JSON object.
    """
    #first checks if the data is a dictionary then checks every key:value's value to see if that is a dictionary
    if isinstance(data, dict):
        return {k: sanitise_json(v) for k, v in data.items()}
    #then checks for lists
    elif isinstance(data, list):
        return [sanitise_json(v) for v in data]
    #finally checks for strings and then santises the data
    elif isinstance(data, str):
        # Remove control characters (non-printable that can corrupt files etc) and escape characters when it sees one or more of them
        data = re.sub(r'[\x00-\x1F\x7F\r\n\t]+', '', data)
        # Escape HTML entities
        return html.escape(data)
    else:
        return data
    
#type out expected keys and value types, eg. expected_keys = {"username": str, "age": int}

def validate_json_structure(safe_data, expected_keys):
    valid = True
    for key, expected_type in expected_keys.items():
        if key not in safe_data:
            print(f"Error: Missing required key")
            valid = False
        elif not isinstance(safe_data[key], expected_type):
            print(f"Error: Invalid key types")
            valid = False
    return valid

def apitest(api_url):
    '''
    Tests an api url to see if a connection can be established and if the response is valid json.
    Then it sanitises the json data for anything malicious before printing the data as a python dict
    '''
    try:
        response = requests.get(api_url, timeout=10)
        if response.status_code == 200:
            try:
                data = response.json()
            except ValueError:
                print("Error: Reposnse is not valid JSON")
                return
            
            safe_data = sanitise_json(data)
            if validate_json_structure(safe_data, expected_keys): 
                print(json.dumps(safe_data, indent = 4))
            else:
                print("Invalid or missing key")
        else:
            print("Error: Response code =", response.status_code)
    except requests.RequestException as e:
        print("Error: Could not connect to API", e)
        return


#to use the code
expected_keys = {"login": str, "updated_at": str}
api_endpoint = "https://api.github.com/users/octocat"
#or
expected_keys = {"craft": str, "name": str}
api_endpoint = "http://api.open-notify.org/astros.json"

apitest(api_endpoint)