import json
import os

data = []
for file in os.listdir('/home/vaughan/rsz/'):
	if file.startswith('xid') and file.endswith('.json'):
		with open(file) as json_file:
			datastore = json.load(json_file)
			data.append(datastore)

for key in data[0].keys():
	print(data[0][key])
