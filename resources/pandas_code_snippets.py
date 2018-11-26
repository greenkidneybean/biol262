# code snippets

# import pandas
import pandas as pd

# create dataframe
df = pd.DataFrame('data/ILINet.csv')

# functions to explore data
df.shape # returns the number of rows and columns
df.info() # returns information regarding column types and null values
df.describe() # returns general descriptive statistics for numeric columns
df.columns # returns the string labels for each column
df.head() # returns the top 5 rows of the dataframe
df.tail() # returns the last 5 rows of the dataframe
df.sample(10) # returns a random sampling of rows from the dataframe

# wrangling functions
df.isnull().sum() # combining 2 functions to get the count of null values in each column
df.replace() # replace a value in a column, or the name of a column
df['column_name'].astype('int') # changes the type of the column ('int', 'float','object')
df['column_name'].unique() # returns list of all unique values in column
df['column_name'].nunique() # returns number of unique values in column
df.pivot() # return a reshaped dataframe organized by provided index/column values

# saving a dataframe:
df.to_csv('file_name.csv') # save your dataframe as a comma separated file

# visualizing dataframes
# helpful packages and settings
import matplotlib.pyplot as plt
import seaborn as sns
plt.style.use('seaborn')
%matplotlib inline # this is specific to jupyter notebooks

# simple plots using Pandas
df.plot()
df.hist()
df.boxplot()

# plots using Seaborn
plt = sns.swarmplot()
fig = fig.get_figure()
fig.save_fig('name_of_figure')
