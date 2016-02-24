
# coding: utf-8

# In[63]:

import csv; import re
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import cluster, decomposition


# In[145]:

data = pd.read_csv("../data/all_condensed_v3_nk.csv", error_bad_lines=False, encoding='cp1252')


# In[29]:

#data.loc[data['Country'] == "Cote d'Ivoire"]
#data.replace({'C™te dÕIvoire': "Cote de'ivoire"}, regex=True)
#data['Country'].replace('C(.{,8}voire)', "Cote de'ivoire", inplace=True, regex=True)
#data.loc[data['Country'] == "Cote de'ivoire"]


# In[146]:

# Clean up the dataset
data['Per_TBwithHIV_tested'].replace("&lt;0.1", 0.09, inplace=True)
data.replace("<0.1", 0.09, inplace=True)
data.Prevalence_male_smoking = data.Prevalence_male_smoking.str.replace(' ', '')
data.Prevalence_male_smoking = data.Prevalence_male_smoking.str.replace("\[.*?\]", '')
data.Prevalence_female_smoking = data.Prevalence_female_smoking.str.replace(' ', '')
data.Prevalence_female_smoking = data.Prevalence_female_smoking.str.replace("\[.*?\]", '')


# In[155]:

data = data.convert_objects(convert_numeric=True); data.dtypes
data.replace(r'\s+', np.nan, regex=True)
data.to_csv('../data/all_condensed_v4.csv')


# In[48]:

# Transform dataframe into matrix with n countries and m measurements
len(data.columns)
raw = data.ix[:,2:26]
raw_wona = raw.dropna()
raw_mat = raw.as_matrix()
raw_mat_trans = raw_mat.T; raw_mat_trans.shape


# In[45]:

# Now for PCA ---- work in progress ----
pca = decomposition.PCA()
pca.fit(raw_mat_trans)


# In[ ]:



