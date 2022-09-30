#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @Time   : 2022/9/23
# @Author : Dai Yulong

import os
import pandas as pd

os.chdir(r"E:\projectResearch\02磷酸化蛋白组学流程\parameterResult\1.Info")
df = pd.read_excel('ProbabilitySite.xlsx')

samples = ["PTM.CollapseKey", "window"]
print(samples)
dfsub = df.loc[:, samples]
dfsub.to_csv("Motif_All.txt",sep="\t",index=False)

dfID = dfsub.iloc[0:, -1:]
dfID.to_csv("All_MotifInput.txt",sep="\t",index=False,header=False)