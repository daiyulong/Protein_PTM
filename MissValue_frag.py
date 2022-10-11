#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @Time   : 2022/10/11
# @Author : Dai Yulong

import os
import pandas as pd
import numpy as np
from sklearn.impute import KNNImputer
"""
针对是MSFragger搜库的修饰位点表结果文件
Protein Position\tGene\tProtein\tPeptide\twindow\tBest.Localized.Probability
"""
os.chdir(r"E:\projects_2022\BP20220597-1")
df = pd.read_excel('PTM_site.xlsx')
gdf = pd.read_excel('parameter.xlsx',sheet_name='Sheet1')

samples = ["Protein Position"] + gdf.loc[:, 'Sample'].tolist()
print(samples)
dfsub = df.loc[:, samples]


### 将样本分组文件转换为dict
d_group = {}
for i in range(len(gdf)):
    samp = gdf.iloc[i,0]
    group = gdf.iloc[i,1]
    if group not in d_group:
        d_group[group] = [samp]
    else:
        d_group[group].append(samp)
print(d_group)
df_data_copy = dfsub
df_data_copy = df_data_copy.replace(0, np.NaN)   # 将intensity中的0值替换成NA值

"""缺失值过滤"""
### 所有组内样本缺失值大于50%的物质滤除
def Confidence(df=None, d_group=None, confidence=0.5,how="outer",):
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
    df_filter.to_csv("missValue_filter.csv",index=0)
    return df_


df2 = Confidence(df=df_data_copy,d_group=d_group,confidence=0.5,how='outer')
df2.to_excel("missValue_clean.xlsx",index=None)


## KNN 方法进行缺失值填充
imputer = KNNImputer(n_neighbors=10)
df2.iloc[:,1:] = imputer.fit_transform(df2.iloc[:,1:])

print(df2.head(5).to_string())
df2.to_excel("missValue_imputation.xlsx",index=None)
