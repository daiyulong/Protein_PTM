import os
import pandas as pd
import numpy as np
from sklearn.impute import KNNImputer


os.chdir(r"E:\projects_2022\BP20220448\infile\codeTest")
df = pd.read_excel('AllSites.xlsx')
gdf = pd.read_excel('parameter.xlsx',sheet_name='Sheet1')
comdf = pd.read_excel('parameter.xlsx',sheet_name='Sheet2')

samples =  ["PTM.CollapseKey"] + gdf.loc[:,'Sample'].tolist()
dfsub = df.loc[:,samples]
dfsub.set_index('PTM.CollapseKey')
#print(dfsub.head().to_string())

gdf.to_csv("groupFile.txt",sep="\t",index=False)


# 将样本分组文件转换为字典形式
d_group = {}
for i in range(len(gdf)):
    samp = gdf.iloc[i,0]
    group = gdf.iloc[i,1]
    if group not in d_group:
        d_group[group] = [samp]
    else:
        d_group[group].append(samp)
print(d_group)
df_data_copy = dfsub.copy(deep=True)
df_data_copy = df_data_copy.replace('Filtered', np.NaN)


def Confidence2(df=None, d_group=None, confidence=0.5,how="outer",):
    """
    函数接收一个数据框和分组，实现自定义分组过滤,
    当confidence=0.5 时， how='outer',表示但凡有一组满足confidence时保留，所有组均不满足confidence时过滤掉。
    当confidence=0.5 时， how='inner',表示需所有组满足confidence时保留，但凡有一组满足不满足confidence时过滤掉
    :param df: 需要处理的数据框
    :param d_group: 数据框中的样本分组设置
    """
    d_set = set()
    l_list = []
    for i in d_group:
        dx = df.loc[:, d_group[i]]
        n = dx.shape[1] * confidence
        print(n)
        ind = dx[dx.isna().sum(axis=1) < dx.shape[1] - n].index

        [d_set.add(i) for i in ind]
        l_list.append(list(ind))

    from functools import reduce
    intersection = sorted(list(reduce(lambda x, y: set(x) & set(y), l_list)))
    union = sorted(list(d_set))

    if how == "inner":
        df_ = df.loc[intersection,:]
    elif how == "outer":
        df_ = df.loc[union,:]

    df_filter = df.loc[[x for x in df.index if x not in df_.index],:]
    df_filter.to_csv("filter.csv",index=0)
    return df_


df2 = Confidence2(df=df_data_copy,d_group=d_group,confidence=0.5,how='outer')
df2.to_excel("缺失值过滤后数据.xlsx",index=None)


# ???? 参数控制：KNN填充，最小值填充，均值填充等等
## KNN method imputer
imputer = KNNImputer(n_neighbors=10)
df2.iloc[:,1:] = imputer.fit_transform(df2.iloc[:,1:])

print(df2.head(5).to_string())
df2.to_excel("KNN填充后数据.xlsx",index=None)



