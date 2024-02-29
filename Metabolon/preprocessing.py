# class to preprocess metabolon data

import pandas as pd
import numpy as np
import os
import re
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

class MetabolonDataset:
    '''
    Class to load and QC metabolon data
    '''
    def __init__(self, file_path, id, node_name):
        self.file_path = file_path
        self.raw_data = None
        self.compound_mappers = None
        self.processed_data = None
        self.metadata = None
        self.id = id
        self.node_name = node_name

        self.read_data(file_path)
        
    def read_data(self, file_path):
        """
        Read in the metabolon data
        """
        # read in the data
        data = pd.read_excel(file_path + '/peaktable.xlsx')
        self.raw_data = data

        metadata = pd.read_csv(file_path + '/s_' + self.id + '.txt', sep = '\t')

        self.metadata = metadata

        first_col_name = data.columns[data.columns.str.contains('Unnamed') == False][0]
        first_col_index = data.columns.get_loc(first_col_name)

        self.compound_mappers = data.iloc[:, 0:first_col_index]
        
        data_filt = data.iloc[3:, :]
        data_filt.columns = data.iloc[2]
        data_filt = data_filt.iloc[:, first_col_index:]
        print(data_filt)


        return data
    
    def preprocess_data(self):
        pass

pmh_data = MetabolonDataset(file_path = 'MTBLS136', id = 'MTBLS136', node_name = 'PMH')

print(pmh_data)