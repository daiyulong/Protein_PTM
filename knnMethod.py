import numpy as np
import pandas as pd
import pandas_profiling as pp
import matplotlib.pyplot as plt
import seaborn as sns
import os
sns.set(context="notebook", style="darkgrid")

from sklearn.impute import KNNImputer

os.chdir(r"E:\projects_2022\20220809")
df = pd.read_excel('pep-new.xlsx',index_col=0)

df.columns = ['[1] PLC-NC1.PTM.Quantity','[2] PLC-NC2.PTM.Quantity','[3] PLC-NC3.PTM.Quantity',
                    '[4] PLC-OE1.PTM.Quantity','[5] PLC-OE2.PTM.Quantity','[6] PLC-OE3.PTM.Quantity']
headline = df.head()
print(headline.to_string())
#查看每个样本中的缺失值个数
df_data_copy = df.copy(deep=True)
df_data_copy[['[1] PLC-NC1.PTM.Quantity','[2] PLC-NC2.PTM.Quantity','[3] PLC-NC3.PTM.Quantity',
                    '[4] PLC-OE1.PTM.Quantity','[5] PLC-OE2.PTM.Quantity','[6] PLC-OE3.PTM.Quantity']] = df_data_copy[['[1] PLC-NC1.PTM.Quantity','[2] PLC-NC2.PTM.Quantity','[3] PLC-NC3.PTM.Quantity',
                    '[4] PLC-OE1.PTM.Quantity','[5] PLC-OE2.PTM.Quantity','[6] PLC-OE3.PTM.Quantity']].replace('Filtered', np.NaN)
print(df_data_copy.isnull().sum())
#查看样本中缺失值的位置
null_index = df_data_copy.loc[df_data_copy['[1] PLC-NC1.PTM.Quantity'].isnull(), :].index
print(null_index)
#
imputer = KNNImputer(n_neighbors=10)
df_data_copy[['[1] PLC-NC1.PTM.Quantity', '[2] PLC-NC2.PTM.Quantity', '[3] PLC-NC3.PTM.Quantity',
              '[4] PLC-OE1.PTM.Quantity','[5] PLC-OE2.PTM.Quantity','[6] PLC-OE3.PTM.Quantity']] = imputer.fit_transform(df_data_copy[['[1] PLC-NC1.PTM.Quantity', '[2] PLC-NC2.PTM.Quantity', '[3] PLC-NC3.PTM.Quantity',
                                                                                                                                             '[4] PLC-OE1.PTM.Quantity','[5] PLC-OE2.PTM.Quantity','[6] PLC-OE3.PTM.Quantity']])
print(df_data_copy.isnull().sum())
headline2 = df_data_copy.head()
print(headline2.to_string())

df_data_copy.to_excel("KnnImpute.xlsx")
