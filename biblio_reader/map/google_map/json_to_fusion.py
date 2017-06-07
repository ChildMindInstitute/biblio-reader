import manager as mg, json, os, pandas as pd

affiliations = json.load(mg.get_file('affiliations.json', os.curdir))

data = mg.get_data()

fusion = pd.Series({latlong: ', '.join([data.loc[i, 'Title'] for i in affiliations[latlong]["papers"]])
                    for latlong in affiliations}, name='Titles')

fusion.to_csv(path='fusion.csv', index_label='Location', header=True)

