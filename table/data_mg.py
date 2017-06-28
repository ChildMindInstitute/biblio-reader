# Used to isolate specific items from the Pandas data structure
# Can be isolated based on publication type, whether it uses the dataset, the affiliation, author, etc.

import sys
sys.path.insert(0, "/Users/jake.son/PycharmProjects/Biblio_Reader")
import manager as mg
data = mg.get_data()

# Potentially create a function for this

data = data[data['Data Use'] == 'Y'][data['Journal Category'] == 'Journal']
table_data = data[['Title', 'Authors', 'Year', 'URL']].sort_values(['Year', 'Title'])
table_data.to_csv(path_or_buf='Table_data.csv', index=False)