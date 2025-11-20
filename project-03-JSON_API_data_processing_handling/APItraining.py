
#Status codes
#200 OK: Request successful, data returned.
#201 Created: New resource created.
#204 No Content: Success but no data returned.
#400 Bad Request: Invalid request.
#401 Unauthorized: Missing or invalid API key.
#500 Internal Server Error: Server encountered an error.


#this test prints the deserielised data from an api as a python dict to inspect before parsing

import requests
import json

def apitest():
    url = "https://api.github.com/users/octocat"
    response = requests.get(url)
    data = response.json()
    print(json.dumps(data, indent = 4))

#same test but api url can be changed by the user

import requests
import json

def apitest(api_url):
    response = requests.get(api_url)
    data = response.json()
    print(json.dumps(data, indent = 4))

api_endpoint = "https://api.github.com/users/octocat"
#or
api_endpoint = "http://api.open-notify.org/astros.json"

apitest(api_endpoint)


#get location data for any user via the github user api

def get_location_data(api_url):
    response = requests.get(api_url)
    if response.status_code == 200:
        data = response.json()
        last_updated = data["updated_at"]
        location = data["location"]
        return location, last_updated
    else:
        return None

api_endpoint = "https://api.github.com/users/octocat"
requested_info = get_location_data(api_endpoint)
symbol = "Requested User Information:"
if requested_info is not None:
    print(f"{symbol}: {requested_info}")
else:
    print("Failed to retreive data")