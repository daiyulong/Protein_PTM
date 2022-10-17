#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @Time   : 2022/9/19
# @Author : Dai Yulong

import argparse
import os
import sys
import shutil
import platform
import pandas as pd
import numpy as np
import subprocess
from sklearn.impute import KNNImputer
import wfmo


"""MissValue Analysis"""
def miss_value_handling(input1,para,pn,isft):
    of = pn + "/1.Info"
    if not os.path.exists(of):
        os.mkdir(of)
    if os.path.exists(input1):
        df = pd.read_excel(input1)
        gdf = pd.read_excel(para, sheet_name='Sheet1')
        if isft == "s":
            samples = ["PTM.CollapseKey"] + gdf.loc[:, 'Sample'].tolist()
            print(samples)
            dfsub = df.loc[:, samples]
            dfsub.set_index('PTM.CollapseKey')

            # 通过parameter.xlsx获得样本及分组信息
            d_group = {}
            for i in range(len(gdf)):
                samp = gdf.iloc[i, 0]
                group = gdf.iloc[i, 1]
                if group not in d_group:
                    d_group[group] = [samp]
                else:
                    d_group[group].append(samp)
            print(d_group)
            df_data_copy = dfsub
            df_data_copy = df_data_copy.replace('Filtered', np.NaN)
        elif isft == "m":
            samples = ["Protein Position"] + gdf.loc[:, 'Sample'].tolist()
            print(samples)
            dfsub = df.loc[:, samples]

            ### 将样本分组文件转换为dict
            d_group = {}
            for i in range(len(gdf)):
                samp = gdf.iloc[i, 0]
                group = gdf.iloc[i, 1]
                if group not in d_group:
                    d_group[group] = [samp]
                else:
                    d_group[group].append(samp)
            print(d_group)
            df_data_copy = dfsub
            df_data_copy = df_data_copy.replace(0, np.NaN)  # 将intensity中的0值替换成NA值

        def Confidence(df=None, d_group=None, confidence=0.5, how="outer", ):
            """函数接收一个数据框和分组，实现自定义分组过滤"""
            d_set = set()
            l_list = []
            for i in d_group:
                dx = df.loc[:, d_group[i]]
                n = dx.shape[1] * confidence
                #print(n)
                ind = dx[dx.isna().sum(axis=1) < dx.shape[1] - n].index
                [d_set.add(i) for i in ind]
                l_list.append(list(ind))

            from functools import reduce
            intersection = sorted(list(reduce(lambda x, y: set(x) & set(y), l_list)))
            union = sorted(list(d_set))

            if how == "inner":
                df_ = df.loc[intersection, :]
            elif how == "outer":
                df_ = df.loc[union, :]

            df_filter = df.loc[[x for x in df.index if x not in df_.index], :]
            df_filter.to_csv(of + "/missValue_filter.csv", index=0)
            return df_

        df2 = Confidence(df=df_data_copy, d_group=d_group, confidence=0.5, how='outer')
        df2.to_excel(of + "/missValue_clean.xlsx", index=None)

        ## KNN method imputer
        imputer = KNNImputer(n_neighbors=10)
        df2.iloc[:, 1:] = imputer.fit_transform(df2.iloc[:, 1:])

        print(df2.head(5).to_string())
        df2.to_excel(of + "/missValue_imputation.xlsx", index=None)


"""Sample Analysis"""
def run_sample_analysis(groupFile, scaleMethod, col, pn, bin, olog, delmin, exp='DIA'):
    of = pn + "/2.SampleAnalysis"
    if not os.path.exists(of):
        os.mkdir(of)
    if os.path.exists(groupFile):

        # 探索性绘图
        print("### Run matrix drawing ###")
        olog.write(
            "Rscript {0}/matrixDraw.R {1} {2} {3} {4} Intensity log2 {5}\n".format(
                bin, pn + "/1.Info/matrix.txt", groupFile, of, col, delmin))
        os.system(
            "Rscript {0}/matrixDraw.R {1} {2} {3} {4} Intensity log2 {5}".format(
                bin, pn + "/1.Info/matrix.txt", groupFile, of, col, delmin))

        olog.write("Rscript {}/svdiPCA2.R {} {} {} uv\n".format(bin, pn + "/1.Info/matrix.txt", groupFile, of))
        os.system("Rscript {}/svdiPCA2.R {} {} {} uv".format(bin, pn + "/1.Info/matrix.txt", groupFile, of))

"""获得比较组信息"""
def get_comparison(infile):
    ifile = open(infile)
    L = dict()
    for line in ifile.readlines()[1:]:
        row = line.rstrip().split('\t')
        row[0] = row[0].lower()
        if row[0] == 'ttest' or row[0] == 'pairedttest':
            L["{}.vs.{}".format(row[1], row[2])] = row[0]
    return (L)

"""通过group文件的第二列获取group数"""
def get_group_numb(groupfile):
    ifile = open(groupfile)
    L = []
    for line in ifile.readlines()[1:]:
        row = line.rstrip().split('\t')
        L.append(row[1])
    return (len(set(L)))

"""获取所有可信位点的motif序列"""
def get_all_motif_seq(probfile, pn, isft):
    of = pn + "/MotifAnalysis"
    if not os.path.exists(of):
        os.mkdir(of)
    if os.path.exists(probfile):
        df = pd.read_excel(probfile)
        if isft == "s":
            samples = ["PTM.CollapseKey", "window"]
            print(samples)
            dfsub = df.loc[:, samples]
            dfsub.to_csv(of + "/Motif_All.txt", sep="\t", index=False)
            dfID = dfsub.iloc[0:, -1:]
            dfID.to_csv(of + "/All_MotifInput.txt", sep="\t", index=False, header=False)
        elif isft == "m":
            samples = ["Protein.Position", "window"]
            print(samples)
            dfsub = df.loc[:, samples]
            dfsub.to_csv(of + "/Motif_All.txt", sep="\t", index=False)
            dfID = dfsub.iloc[0:, -1:]
            dfID.to_csv(of + "/All_MotifInput.txt", sep="\t", index=False, header=False)

"""全部蛋白的分析"""
def run_all_protein_analysis(bin, pn, org, celllocation):

    olog.write("Rscript {}/FuncAnal4ID.R -i {} -o {} -s {} -p {} -c 4\n".format(bin, pn + '/input.txt', pn + '/AllProtein', org, bin))
    os.system("Rscript {}/FuncAnal4ID.R -i {} -o {} -s {} -p {} -c 4".format(bin, pn + '/input.txt', pn + '/AllProtein', org, bin))

    ########细胞定位的分析的绘图
    if os.path.exists(p.pn + '/AllProtein/GOEnrich_Species/CC.txt'):
        olog.write("Rscript {}/dCellLoc.R {} {} {}\n".format(bin, pn + '/AllProtein/GOEnrich_Species/CC.txt', pn + '/AllProtein/GO_subcellular_localization', celllocation))
        os.system("Rscript {}/dCellLoc.R {} {} {}".format(bin, pn + '/AllProtein/GOEnrich_Species/CC.txt', pn + '/AllProtein/GO_subcellular_localization', celllocation))

    #####增加KEGG分类信息
    olog.write("{}/addURLClassKEGG {}/AllProtein/KEGG/KEGG.enriched.txt {}/AllProtein/KEGG/KEGG.enriched.xls\n".format(bin,pn, pn))
    os.system("{}/addURLClassKEGG {}/AllProtein/KEGG/KEGG.enriched.txt {}/AllProtein/KEGG/KEGG.enriched.xls".format(bin, pn,pn))

    olog.write("Rscript {}/txt2xlsx.R {}/AllProtein/KEGG/KEGG.enriched.xls KEGG {}/AllProtein/KEGG/KEGG.enriched.xlsx\n".format(bin, pn, pn))
    os.system("Rscript {}/txt2xlsx.R {}/AllProtein/KEGG/KEGG.enriched.xls KEGG {}/AllProtein/KEGG/KEGG.enriched.xlsx".format(bin, pn, pn))

    #####所有蛋白结果的绘图
    olog.write("Rscript {}/PAAdraw.R {}\n".format(bin, pn + '/AllProtein'))
    os.system("Rscript {}/PAAdraw.R {}".format(bin, pn + '/AllProtein'))

    ############挑选出GO的对应表
    olog.write("Rscript {}/selectGO.R {} 4 {} {}\n".format(bin, pn + '/input.txt', org, pn + "/AllProtein/GO"))
    os.system("Rscript {}/selectGO.R {} 4 {} {}".format(bin, pn + '/input.txt', org, pn + "/AllProtein/GO"))

    ###########生成Summary文件
    olog.write("{}/makeSummary -i {}\n".format(bin, pn))
    os.system("{}/makeSummary -i {}".format(bin, pn))

    olog.write("Rscript {}/tsv2xlsx.R {} {} {}\n".format(bin, pn + '/protein_summary.txt', 'summary', pn + '/summary.xlsx'))
    os.system("Rscript {}/tsv2xlsx.R {} {} {}".format(bin, pn + '/protein_summary.txt', 'summary', pn + '/summary.xlsx'))


## ????
"""多组之间比较的方差分析"""
def run_anova_analysis(bin, pn, org, celllocation, input_file, input_group_file, grouporder):
    olog.write("Rscript {}/extractTable.R {} {}\n".format(bin, input_file, pn + '/1.Info/matrix.txt'))
    os.system("Rscript {}/extractTable.R {} {}".format(bin, input_file, pn + '/1.Info/matrix.txt'))

    olog.write("Rscript {}/kmeans.R {} {} {} {}\n".format(bin, input_group_file, pn + '/1.Info/matrix.txt', pn + '/kmeans', grouporder))
    os.system("Rscript {}/kmeans.R {} {} {} {}".format(bin, input_group_file, pn + '/1.Info/matrix.txt', pn + '/kmeans', grouporder))

    # 挑选做比较绘图的文件
    olog.write("Rscript {}/selectFileFromTotal.R {} {} {}\n".format(bin, pn + '/kmeans/ANOVA_result.xlsx', pn + '/input.txt', pn + '/ANOVA.sig.txt'))
    os.system("Rscript {}/selectFileFromTotal.R {} {} {}".format(bin, pn + '/kmeans/ANOVA_result.xlsx', pn + '/input.txt', pn + '/ANOVA.sig.txt'))

    # 比较绘图
    olog.write("Rscript {}/dDEP4multi.R {} {} {}\n".format(bin, pn + '/ANOVA.sig.txt', input_group_file, pn + '/Multi_group_expression'))
    os.system("Rscript {}/dDEP4multi.R {} {} {}".format(bin, pn + '/ANOVA.sig.txt', input_group_file, pn + '/Multi_group_expression'))

    #热图绘制
    olog.write("Rscript {}/drawHeatmap.R {} {} {} multi\n".format(bin, pn + '/ANOVA.sig.txt',input_group_file, pn))
    os.system("Rscript {}/drawHeatmap.R {} {} {} multi".format(bin, pn + '/ANOVA.sig.txt',input_group_file, pn))

    # 显著的蛋白的功能分析和绘图
    olog.write("Rscript {}/FuncAnal4ID.R -i {} -o {} -s {} -p {} -b {} -d 4 -c 4\n".format(bin, pn + '/ANOVA.sig.txt', pn, org, bin, input_file))
    os.system("Rscript {}/FuncAnal4ID.R -i {} -o {} -s {} -p {} -b {} -d 4 -c 4".format(bin, pn + '/ANOVA.sig.txt', pn, org, bin, input_file))

    ########细胞定位的分析的绘图
    if os.path.exists(p.pn + '/AllProtein/GOEnrich_Species/CC.txt'):
        olog.write("Rscript {}/dCellLoc.R {} {} {}\n".format(bin, pn + '/GOEnrich_Species/CC.txt',pn + '/GO_subcellular_localization', celllocation))
        os.system("Rscript {}/dCellLoc.R {} {} {}".format(bin, pn + '/GOEnrich_Species/CC.txt',pn + '/GO_subcellular_localization', celllocation))

    #####增加KEGG分类信息
    olog.write("{}/addURLClassKEGG {}/KEGG/KEGG.enriched.txt {}/KEGG/KEGG.enriched.xls\n".format(bin, pn, pn))
    os.system("{}/addURLClassKEGG {}/KEGG/KEGG.enriched.txt {}/KEGG/KEGG.enriched.xls".format(bin, pn, pn))

    olog.write("Rscript {}/txt2xlsx.R {}/KEGG/KEGG.enriched.xls KEGG {}/KEGG/KEGG.enriched.xlsx\n".format(bin, pn, pn))
    os.system("Rscript {}/txt2xlsx.R {}/KEGG/KEGG.enriched.xls KEGG {}/KEGG/KEGG.enriched.xlsx".format(bin, pn, pn))

    #####结果的绘图
    olog.write("Rscript {}/PAAdraw.R {}\n".format(bin, pn))
    os.system("Rscript {}/PAAdraw.R {}".format(bin, pn))

    if os.path.exists(pn + '/ANOVA.sig.txt'):
        ih = open(pn + '/ANOVA.sig.txt')
        lines = ih.readlines()
        return len(lines[1:])
    else:
        return 0

def get_python():
    if platform.system()=='Windows':
        return ("python")
    else:
        return("python3")

"""差异Motif分析"""
def get_DEP_motif(L,pn):
    for comp in L:
        print("\nRun " + comp)
        infile = pn + "/MotifAnalysis/Motif_All.txt"
        DEPfile = pn + '/' + comp + '/DESelection/DEP.txt'
        if os.path.exists(infile) and os.path.exists(DEPfile):
            '''提取DEP Motif分析的输入序列 '''
            olog.write("Rscript {}/getDEPMotif.R {} {}\n".format(bin, comp, pn))
            os.system("Rscript {}/getDEPMotif.R {} {}".format(bin, comp, pn))

            # DEP motif分析
            of = os.path.join(os.getcwd(), p.pn + '/DEPMotif/' + comp)
            print(of)
            if not os.path.exists(of):
                os.mkdir(of)

            if not os.path.exists(of + '/All'):
                os.mkdir(of + '/All')
            if os.path.exists(of + '/DEPMotif_input.txt'):
                codeContent = 'momo motifx -oc /meme/All --verbosity 1 --width 13 --eliminate-repeats 13 --min-occurrences 5'
                subprocess.run(
                    "echo \"welovebp1188\" |sudo -S docker run -v {workpath}:/meme --user `id -u`:`id -g` {docker_name} {codeC}  /meme/DEPMotif_input.txt"
                    .format(workpath=of, docker_name='memesuite/memesuite:5.3.3', codeC=codeContent), shell=True,
                    check=True)
            else:
                print("未检测到DEPMotif_input.txt输入文件！")
                sys.exit()

            # Up motif分析
            if not os.path.exists(of + '/Up'):
                os.mkdir(of + '/Up')
            if os.path.exists(of + '/UpMotif_input.txt'):
                codeContent = 'momo motifx -oc /meme/Up --verbosity 1 --width 13 --eliminate-repeats 13 --min-occurrences 5'
                subprocess.run(
                    "echo \"welovebp1188\" |sudo -S docker run -v {workpath}:/meme --user `id -u`:`id -g` {docker_name} {codeC}  /meme/UpMotif_input.txt"
                    .format(workpath=of, docker_name='memesuite/memesuite:5.3.3', codeC=codeContent), shell=True,
                    check=True)
            else:
                print("未检测到UpMotif_input.txt输入文件！")

            # Down motif分析
            if not os.path.exists(of + '/Down'):
                os.mkdir(of + '/Down')
            if os.path.exists(of + '/DownMotif_input.txt'):
                codeContent = 'momo motifx -oc /meme/Down --verbosity 1 --width 13 --eliminate-repeats 13 --min-occurrences 5'
                subprocess.run(
                    "echo \"welovebp1188\" |sudo -S docker run -v {workpath}:/meme --user `id -u`:`id -g` {docker_name} {codeC}  /meme/DownMotif_input.txt"
                        .format(workpath=of, docker_name='memesuite/memesuite:5.3.3', codeC=codeContent),
                    shell=True,
                    check=True)
            else:
                print("未检测到DownMotif_input.txt输入文件！")

"""差异蛋白list生成"""
def get_DEP_protein_id(L,pn,isft):
    for comp in L:
        print("\nRun " + comp + "\n")
        infile = pn + "/1.Info/ModifiedProtein.xlsx"
        DEPfile = pn + "/input.txt"
        if os.path.exists(infile) and os.path.exists(DEPfile):
            olog.write("Rscript {}/getDEPprotein.R {} {} {}\n".format(bin, comp, pn, isft))
            os.system("Rscript {}/getDEPprotein.R {} {} {}".format(bin, comp, pn, isft))

"""差异位点的筛选"""
def run_dep_selection(L, bin, pn, olog, fc=2, pvalue=0.05):
    for comp in L:
        infile = pn + '/' + comp + "/input.txt"
        groupfile = pn + '/' + comp + '/group.txt'
        if os.path.exists(infile) and os.path.exists(groupfile):
            DEP = pn + '/' + comp + '/DESelection'
            # 提取数据
            olog.write("Rscript {}/extractTablePTM.R {} {}\n".format(bin, infile, pn + '/' + comp + '/matrix.txt'))
            os.system("Rscript {}/extractTablePTM.R {} {}".format(bin, infile, pn + '/' + comp + '/matrix.txt'))
            # 比较组之间绘图不进行

            # 差异筛选
            if L[comp] == 'ttest':
                olog.write("Rscript {}/Ttest.R {}/matrix.txt {} {} none {} {}\n".format(bin, pn + '/' + comp, DEP, groupfile, fc, pvalue))
                os.system("Rscript {}/Ttest.R {}/matrix.txt {} {} none {} {}".format(bin, pn + '/' + comp, DEP, groupfile, fc, pvalue))
            elif L[comp] == 'pairedttest':
                olog.write("Rscript {}/Ttest.R {}/matrix.txt {} {} none {} {} TRUE\n".format(bin, pn + '/' + comp, DEP, groupfile, fc, pvalue))
                os.system("Rscript {}/Ttest.R {}/matrix.txt {} {} none {} {} TRUE".format(bin, pn + '/' + comp, DEP, groupfile, fc, pvalue))

            # 合并数据
            olog.write("Rscript {}/mergeTablePTM.R {}/DEP.txt {} {}/DEP.xls\n".format(bin, DEP, infile, DEP))
            os.system("Rscript {}/mergeTablePTM.R {}/DEP.txt {} {}/DEP.xls".format(bin, DEP, infile, DEP))

            olog.write("Rscript {}/mergeTablePTM.R {}/ttest.txt {} {}/AllTest.xls\n".format(bin, DEP, infile, DEP))
            os.system("Rscript {}/mergeTablePTM.R {}/ttest.txt {} {}/AllTest.xls".format(bin, DEP, infile, DEP))

            olog.write("Rscript {}/mergeTablePTM.R {}/Down.txt {} {}/Down.xls\n".format(bin, DEP, infile, DEP))
            os.system("Rscript {}/mergeTablePTM.R {}/Down.txt {} {}/Down.xls".format(bin, DEP, infile, DEP))

            olog.write("Rscript {}/mergeTablePTM.R {}/Up.txt {} {}/Up.xls\n".format(bin, DEP, infile, DEP))
            os.system("Rscript {}/mergeTablePTM.R {}/Up.txt {} {}/Up.xls".format(bin, DEP, infile, DEP))


            # 差异绘图，火山图，差异蛋白统计图
            olog.write("Rscript {}/dVolcano.R {}/AllTest.xls {} {} {}\n".format(bin, DEP, fc, pvalue, DEP))
            os.system("Rscript {}/dVolcano.R {}/AllTest.xls {} {} {}".format(bin, DEP, fc, pvalue, DEP))

            # 热图绘制
            olog.write("Rscript {}/drawHeatmapPTM.R {}/DEP.xls {} {}\n".format(bin, DEP, groupfile, DEP))
            os.system("Rscript {}/drawHeatmapPTM.R {}/DEP.xls {} {}".format(bin, DEP, groupfile, DEP))

            # 差异比较绘图
            olog.write("Rscript {}/dDEP4PTM.R {}/DEP.xls {} {}\n".format(bin, DEP, groupfile, DEP))
            os.system("Rscript {}/dDEP4PTM.R {}/DEP.xls {} {}".format(bin, DEP, groupfile, DEP))


"""差异位点对应蛋白的功能分析"""
def run_dep_function_analysis(L, bin, pn, org, celllocation, supp, olog):
    for comp in L:
        DEP = pn + '/' + comp + '/DESelection'
        # 提取Summary
        olog.write(
            "Rscript {}/selectTSV.R {}/DEPProtein.xls {}/protein_summary.txt {}/all_summamry.xlsx\n".format(bin, DEP, pn, pn + '/' + comp))
        os.system(
            "Rscript {}/selectTSV.R {}/DEPProtein.xls {}/protein_summary.txt {}/all_summary.xlsx".format(bin, DEP, pn, pn + '/' + comp))

        olog.write(
            "Rscript {}/selectTSV.R {}/UpProtein.xls {}/protein_summary.txt {}/up_summary.xlsx\n".format(bin, DEP, pn, pn + '/' + comp))
        os.system(
            "Rscript {}/selectTSV.R {}/UpProtein.xls {}/protein_summary.txt {}/up_summary.xlsx".format(bin, DEP, pn, pn + '/' + comp))

        olog.write(
            "Rscript {}/selectTSV.R {}/DownProtein.xls {}/protein_summary.txt {}/down_summary.xlsx\n".format(bin, DEP, pn, pn + '/' + comp))
        os.system(
            "Rscript {}/selectTSV.R {}/DownProtein.xls {}/protein_summary.txt {}/down_summary.xlsx".format(bin, DEP, pn, pn + '/' + comp))
        # 功能分析
        print("\nRun " + comp)
        if org in supp:
            of = pn + '/' + comp + "/all"
            ofu = pn + '/' + comp + "/up"
            ofd = pn + '/' + comp + "/down"
            dir_list = [of, ofu, ofd]
            for dir_ in dir_list:
                if not os.path.exists(dir_):
                    os.mkdir(dir_)
            # DEP
            olog.write("Rscript {}/convertID.R {}/DEPProtein.xls {} {} {}\n".format(bin, DEP, of + '/Input', org, p.type))
            os.system("Rscript {}/convertID.R {}/DEPProtein.xls {} {} {}".format(bin, DEP, of + '/Input', org, p.type))
            # Up
            olog.write("Rscript {}/convertID.R {}/UpProtein.xls {} {} {}\n".format(bin, DEP, ofu + '/Input', org, p.type))
            os.system("Rscript {}/convertID.R {}/UpProtein.xls {} {} {}".format(bin, DEP, ofu + '/Input', org, p.type))
            # Down
            olog.write("Rscript {}/convertID.R {}/DownProtein.xls {} {} {}\n".format(bin, DEP, ofd + '/Input', org, p.type))
            os.system("Rscript {}/convertID.R {}/DownProtein.xls {} {} {}".format(bin, DEP, ofd + '/Input', org, p.type))

            if not os.path.exists(of + '/Input/Transformed_input.txt'):
                print("转换DEP ID出错")
                sys.exit()
            if not os.path.exists(ofu + '/Input/Transformed_input.txt'):
                print("转换Up ID出错")
                sys.exit()
            if not os.path.exists(ofd + '/Input/Transformed_input.txt'):
                print("转换Down ID出错")
                sys.exit()

            if (p.ibg == 'none'):
                olog.write("Rscript {}/FuncAnal4ID.R -i {} -s {} -p {} -o {}\n".format(bin, of + '/Input/Transformed_input.txt', org, bin, of))
                os.system("Rscript {}/FuncAnal4ID.R -i {} -s {} -p {} -o {}".format(bin, of + '/Input/Transformed_input.txt', org, bin, of))
                # UP
                olog.write("Rscript {}/FuncAnal4ID.R -i {} -s {} -p {} -o {}\n".format(bin, ofu + '/Input/Transformed_input.txt', org, bin, ofu))
                os.system("Rscript {}/FuncAnal4ID.R -i {} -s {} -p {} -o {}".format(bin, ofu + '/Input/Transformed_input.txt', org, bin, ofu))
                # Down
                olog.write("Rscript {}/FuncAnal4ID.R -i {} -s {} -p {} -o {}\n".format(bin, ofd + '/Input/Transformed_input.txt', org, bin, ofd))
                os.system("Rscript {}/FuncAnal4ID.R -i {} -s {} -p {} -o {}".format(bin, ofd + '/Input/Transformed_input.txt', org, bin, ofd))

            else:
                olog.write("Rscript {}/convertID.R {} {} {} {}\n".format(bin, p.ibg, of + '/BG', org, p.type))
                os.system("Rscript {}/convertID.R {} {} {} {}".format(bin, p.ibg, of + '/BG', org, p.type))
                olog.write("Rscript {}/FuncAnal4ID.R -i {} -s {} -p {} -o {} -b {}\n".format(bin, of + '/Input/Transformed_input.txt',org, bin, of, of + '/BG/Transformed_input.txt'))
                os.system("Rscript {}/FuncAnal4ID.R -i {} -s {} -p {} -o {} -b {}".format(bin, of + '/Input/Transformed_input.txt',org, bin, of, of + '/BG/Transformed_input.txt'))
                #UP
                olog.write("Rscript {}/convertID.R {} {} {} {}\n".format(bin, p.ibg, ofu + '/BG', org, p.type))
                os.system("Rscript {}/convertID.R {} {} {} {}".format(bin, p.ibg, ofu + '/BG', org, p.type))
                olog.write("Rscript {}/FuncAnal4ID.R -i {} -s {} -p {} -o {} -b {}\n".format(bin, ofu + '/Input/Transformed_input.txt',org, bin, ofu, ofu + '/BG/Transformed_input.txt'))
                os.system("Rscript {}/FuncAnal4ID.R -i {} -s {} -p {} -o {} -b {}".format(bin, ofu + '/Input/Transformed_input.txt',org, bin, ofu, ofu + '/BG/Transformed_input.txt'))
                #Down
                olog.write("Rscript {}/convertID.R {} {} {} {}\n".format(bin, p.ibg, ofd + '/BG', org, p.type))
                os.system("Rscript {}/convertID.R {} {} {} {}".format(bin, p.ibg, ofd + '/BG', org, p.type))
                olog.write("Rscript {}/FuncAnal4ID.R -i {} -s {} -p {} -o {} -b {}\n".format(bin, ofd + '/Input/Transformed_input.txt',org, bin, ofd, ofd + '/BG/Transformed_input.txt'))
                os.system("Rscript {}/FuncAnal4ID.R -i {} -s {} -p {} -o {} -b {}".format(bin, ofd + '/Input/Transformed_input.txt',org, bin, ofd, ofd + '/BG/Transformed_input.txt'))

            ########细胞定位的分析的绘图
            #DEP
            if os.path.exists(of + '/GOEnrich_Species/CC.txt'):
                olog.write("Rscript {}/dCellLoc.R {} {} {}\n".format(bin, of + '/GOEnrich_Species/CC.txt', of + '/GO_subcellular_localization', celllocation))
                os.system("Rscript {}/dCellLoc.R {} {} {}".format(bin, of + '/GOEnrich_Species/CC.txt', of + '/GO_subcellular_localization', celllocation))
            #UP
            if os.path.exists(ofu + '/GOEnrich_Species/CC.txt'):
                olog.write("Rscript {}/dCellLoc.R {} {} {}\n".format(bin, ofu + '/GOEnrich_Species/CC.txt', ofu + '/GO_subcellular_localization', celllocation))
                os.system("Rscript {}/dCellLoc.R {} {} {}".format(bin, ofu + '/GOEnrich_Species/CC.txt', ofu + '/GO_subcellular_localization', celllocation))
            #Down
            if os.path.exists(ofd + '/GOEnrich_Species/CC.txt'):
                olog.write("Rscript {}/dCellLoc.R {} {} {}\n".format(bin, ofd + '/GOEnrich_Species/CC.txt', ofd + '/GO_subcellular_localization', celllocation))
                os.system("Rscript {}/dCellLoc.R {} {} {}".format(bin, ofd + '/GOEnrich_Species/CC.txt', ofd + '/GO_subcellular_localization', celllocation))

            #####增加KEGG分类信息
            if p.type != 'accsymbolfc':
                olog.write("{}/addURLClassKEGG {}/KEGG/KEGG.enriched.txt {}/KEGG/KEGG.enriched.xls\n".format(bin, of, of))
                os.system("{}/addURLClassKEGG {}/KEGG/KEGG.enriched.txt {}/KEGG/KEGG.enriched.xls".format(bin, of, of))
                olog.write("Rscript {}/txt2xlsx.R {}/KEGG/KEGG.enriched.xls KEGG {}/KEGG/KEGG.enriched.xlsx\n".format(bin, of, of))
                os.system("Rscript {}/txt2xlsx.R {}/KEGG/KEGG.enriched.xls KEGG {}/KEGG/KEGG.enriched.xlsx".format(bin, of, of))
                #UP
                olog.write("{}/addURLClassKEGG {}/KEGG/KEGG.enriched.txt {}/KEGG/KEGG.enriched.xls\n".format(bin, ofu, ofu))
                os.system("{}/addURLClassKEGG {}/KEGG/KEGG.enriched.txt {}/KEGG/KEGG.enriched.xls".format(bin, ofu, ofu))
                olog.write("Rscript {}/txt2xlsx.R {}/KEGG/KEGG.enriched.xls KEGG {}/KEGG/KEGG.enriched.xlsx\n".format(bin, ofu, ofu))
                os.system("Rscript {}/txt2xlsx.R {}/KEGG/KEGG.enriched.xls KEGG {}/KEGG/KEGG.enriched.xlsx".format(bin, ofu, ofu))
                #Down
                olog.write("{}/addURLClassKEGG {}/KEGG/KEGG.enriched.txt {}/KEGG/KEGG.enriched.xls\n".format(bin, ofd, ofd))
                os.system("{}/addURLClassKEGG {}/KEGG/KEGG.enriched.txt {}/KEGG/KEGG.enriched.xls".format(bin, ofd, ofd))
                olog.write("Rscript {}/txt2xlsx.R {}/KEGG/KEGG.enriched.xls KEGG {}/KEGG/KEGG.enriched.xlsx\n".format(bin, ofd, ofd))
                os.system("Rscript {}/txt2xlsx.R {}/KEGG/KEGG.enriched.xls KEGG {}/KEGG/KEGG.enriched.xlsx".format(bin, ofd, ofd))

            else:
                olog.write("perl {}/addFC4wfsymbol.pl {}\n".format(bin, of))
                os.system("perl {}/addFC4wfsymbol.pl {}".format(bin, of))
                olog.write("Rscript {}/txt2xlsx.R {}/KEGG/KEGG.enriched.addFC.txt KEGG {}/KEGG/KEGG.enriched.xlsx\n".format(bin, of, of))
                os.system("Rscript {}/txt2xlsx.R {}/KEGG/KEGG.enriched.addFC.txt KEGG {}/KEGG/KEGG.enriched.xlsx".format(bin, of, of))
                #UP
                olog.write("perl {}/addFC4wfsymbol.pl {}\n".format(bin, ofu))
                os.system("perl {}/addFC4wfsymbol.pl {}".format(bin, ofu))
                olog.write("Rscript {}/txt2xlsx.R {}/KEGG/KEGG.enriched.addFC.txt KEGG {}/KEGG/KEGG.enriched.xlsx\n".format(bin, ofu, ofu))
                os.system("Rscript {}/txt2xlsx.R {}/KEGG/KEGG.enriched.addFC.txt KEGG {}/KEGG/KEGG.enriched.xlsx".format(bin, ofu, ofu))
                #Down
                olog.write("perl {}/addFC4wfsymbol.pl {}\n".format(bin, ofd))
                os.system("perl {}/addFC4wfsymbol.pl {}".format(bin, ofd))
                olog.write("Rscript {}/txt2xlsx.R {}/KEGG/KEGG.enriched.addFC.txt KEGG {}/KEGG/KEGG.enriched.xlsx\n".format(bin, ofd, ofd))
                os.system("Rscript {}/txt2xlsx.R {}/KEGG/KEGG.enriched.addFC.txt KEGG {}/KEGG/KEGG.enriched.xlsx".format(bin, ofd, ofd))

            #####所有蛋白结果的绘图
            #DEP
            olog.write("Rscript {}/PAAdraw.R {}\n".format(bin, of))
            os.system("Rscript {}/PAAdraw.R {}".format(bin, of))
            #UP
            olog.write("Rscript {}/PAAdraw.R {}\n".format(bin, ofu))
            os.system("Rscript {}/PAAdraw.R {}".format(bin, ofu))
            #Down
            olog.write("Rscript {}/PAAdraw.R {}\n".format(bin, ofd))
            os.system("Rscript {}/PAAdraw.R {}".format(bin, ofd))

            ############挑选出GO的对应表
            #DEP
            olog.write("Rscript {}/selectGO.R {} 1 {} {}\n".format(bin, of + '/Input/Transformed_input.txt', org, of + "/GO"))
            os.system("Rscript {}/selectGO.R {} 1 {} {}".format(bin, of + '/Input/Transformed_input.txt', org, of + "/GO"))
            #UP
            olog.write("Rscript {}/selectGO.R {} 1 {} {}\n".format(bin, ofu + '/Input/Transformed_input.txt', org, ofu + "/GO"))
            os.system("Rscript {}/selectGO.R {} 1 {} {}".format(bin, ofu + '/Input/Transformed_input.txt', org, ofu + "/GO"))
            #Down
            olog.write("Rscript {}/selectGO.R {} 1 {} {}\n".format(bin, ofd + '/Input/Transformed_input.txt', org, ofd + "/GO"))
            os.system("Rscript {}/selectGO.R {} 1 {} {}".format(bin, ofd + '/Input/Transformed_input.txt', org, ofd + "/GO"))

            ###########生成Summary文件
            if p.type == 'symbol':
                olog.write("{}/makeSummary4symbol -i {} -input {}\n".format(bin, of, of + '/Input/Transformed_input.txt'))
                os.system("{}/makeSummary4symbol -i {} -input {}".format(bin, of, of + '/Input/Transformed_input.txt'))
                olog.write("Rscript {}/tsv2xlsx.R {} {} {}\n".format(bin, of + '/protein_summary.txt', 'summary', of + '/summary.xlsx'))
                os.system("Rscript {}/tsv2xlsx.R {} {} {}".format(bin, of + '/protein_summary.txt', 'summary', of + '/summary.xlsx'))
                #UP
                olog.write("{}/makeSummary4symbol -i {} -input {}\n".format(bin, ofu, ofu + '/Input/Transformed_input.txt'))
                os.system("{}/makeSummary4symbol -i {} -input {}".format(bin, ofu, ofu + '/Input/Transformed_input.txt'))
                olog.write("Rscript {}/tsv2xlsx.R {} {} {}\n".format(bin, ofu + '/protein_summary.txt', 'summary', ofu + '/summary.xlsx'))
                os.system("Rscript {}/tsv2xlsx.R {} {} {}".format(bin, ofu + '/protein_summary.txt', 'summary', ofu + '/summary.xlsx'))
                #Down
                olog.write("{}/makeSummary4symbol -i {} -input {}\n".format(bin, ofd, ofd + '/Input/Transformed_input.txt'))
                os.system("{}/makeSummary4symbol -i {} -input {}".format(bin, ofd, ofd + '/Input/Transformed_input.txt'))
                olog.write("Rscript {}/tsv2xlsx.R {} {} {}\n".format(bin, ofd + '/protein_summary.txt', 'summary', ofd + '/summary.xlsx'))
                os.system("Rscript {}/tsv2xlsx.R {} {} {}".format(bin, ofd + '/protein_summary.txt', 'summary', ofd + '/summary.xlsx'))

            ###########PPI search
            taxid = code2taxid[p.org]
            taxid = str(taxid)
            print(taxid)
            if os.path.exists(bin + '/stringdata/' + taxid + '.protein.info.v11.5.txt') and os.path.exists(
                    bin + '/stringdata/' + taxid + '.protein.links.v11.5.txt') and os.path.exists(
                bin + '/stringdata/' + taxid + '.protein.aliases.v11.5.txt'):
                olog.write("{}/ppiSearch -i {} -c {} -t {} -s {} -o {} -pt inter\n".format(bin, of + '/Input/Transformed_input.txt', 3, taxid, 700, of + '/PPI'))
                os.system("{}/ppiSearch -i {} -c {} -t {} -s {} -o {} -pt inter".format(bin, of + '/Input/Transformed_input.txt', 3, taxid, 700, of + '/PPI'))
                olog.write("{} {}/ppi_de-duplicate.py -in {}/PPI.txt -o {}/NonRPPI.txt\n".format(wfmo.get_python(), bin, of + '/PPI', of + '/PPI'))
                os.system("{} {}/ppi_de-duplicate.py -in {}/PPI.txt -o {}/NonRPPI.txt".format(wfmo.get_python(), bin, of + '/PPI', of + '/PPI'))
                #UP
                olog.write("{}/ppiSearch -i {} -c {} -t {} -s {} -o {} -pt inter\n".format(bin, ofu + '/Input/Transformed_input.txt', 3, taxid, 700, ofu + '/PPI'))
                os.system("{}/ppiSearch -i {} -c {} -t {} -s {} -o {} -pt inter".format(bin, ofu + '/Input/Transformed_input.txt', 3, taxid, 700, ofu + '/PPI'))
                olog.write("{} {}/ppi_de-duplicate.py -in {}/PPI.txt -o {}/NonRPPI.txt\n".format(wfmo.get_python(), bin, ofu + '/PPI', ofu + '/PPI'))
                os.system("{} {}/ppi_de-duplicate.py -in {}/PPI.txt -o {}/NonRPPI.txt".format(wfmo.get_python(), bin, ofu + '/PPI', ofu + '/PPI'))
                #Down
                olog.write("{}/ppiSearch -i {} -c {} -t {} -s {} -o {} -pt inter\n".format(bin, ofd + '/Input/Transformed_input.txt', 3, taxid, 700, ofd + '/PPI'))
                os.system("{}/ppiSearch -i {} -c {} -t {} -s {} -o {} -pt inter".format(bin, ofd + '/Input/Transformed_input.txt', 3, taxid, 700, ofd + '/PPI'))
                olog.write("{} {}/ppi_de-duplicate.py -in {}/PPI.txt -o {}/NonRPPI.txt\n".format(wfmo.get_python(), bin, ofd + '/PPI', ofd + '/PPI'))
                os.system("{} {}/ppi_de-duplicate.py -in {}/PPI.txt -o {}/NonRPPI.txt".format(wfmo.get_python(), bin, ofd + '/PPI', ofd + '/PPI'))


                script_dir = os.path.normpath(bin + '/Network')
                model_network_dir = os.path.normpath(script_dir + '/model')
                network_dir = os.path.normpath(of + '/Network')
                unetwork_dir = os.path.normpath(ofu + '/Network')
                dnetwork_dir = os.path.normpath(ofd + '/Network')
                if not os.path.exists(network_dir) or not os.path.exists(unetwork_dir) or not os.path.exists(dnetwork_dir):
                    os.makedirs(network_dir)
                    os.makedirs(unetwork_dir)
                    os.makedirs(dnetwork_dir)

                # 有Fold Change的情况
                if os.path.exists("{}/KEGG/KEGG.enriched.addFC.txt".format(of)) and os.path.exists("{}/KEGG/KEGG.enriched.addFC.txt".format(ofu)) and os.path.exists(
                        "{}/KEGG/KEGG.enriched.addFC.txt".format(ofd)):

                    # DEP
                    os.system("{0} {1}/code_2.py -in {2}/KEGG/KEGG.enriched.addFC.txt -od {3}".format(wfmo.get_python(),script_dir, of,network_dir))
                    olog.write("Rscript {0}/pathway-pathway.R {1}/node_edge/edge_all.txt {1}\n".format(script_dir, network_dir))
                    os.system("Rscript {0}/pathway-pathway.R {1}/node_edge/edge_all.txt {1}".format(script_dir, network_dir))
                    # UP
                    os.system("{0} {1}/code_2.py -in {2}/KEGG/KEGG.enriched.addFC.txt -od {3}".format(wfmo.get_python(),script_dir, ofu,unetwork_dir))
                    olog.write("Rscript {0}/pathway-pathway.R {1}/node_edge/edge_all.txt {1}\n".format(script_dir, unetwork_dir))
                    os.system("Rscript {0}/pathway-pathway.R {1}/node_edge/edge_all.txt {1}".format(script_dir, unetwork_dir))
                    # Down
                    os.system("{0} {1}/code_2.py -in {2}/KEGG/KEGG.enriched.addFC.txt -od {3}".format(wfmo.get_python(),script_dir, ofd,dnetwork_dir))
                    olog.write("Rscript {0}/pathway-pathway.R {1}/node_edge/edge_all.txt {1}\n".format(script_dir, dnetwork_dir))
                    os.system("Rscript {0}/pathway-pathway.R {1}/node_edge/edge_all.txt {1}".format(script_dir, dnetwork_dir))

                    if os.path.exists('{}/node_edge/pathway-pathway_node.txt'.format(network_dir)) and os.path.exists('{}/node_edge/pathway-pathway_node.txt'.format(unetwork_dir)) and os.path.exists(
                            '{}/node_edge/pathway-pathway_node.txt'.format(unetwork_dir)):
                        # DEP
                        olog.write("{0} {1}/pathway-pathway.py -model {2} -net {3} -edge_node {3}\n".format(wfmo.get_python(),script_dir,model_network_dir,network_dir))
                        os.system("{0} {1}/pathway-pathway.py -model {2} -net {3} -edge_node {3}".format(wfmo.get_python(), script_dir, model_network_dir, network_dir))
                        # UP
                        olog.write("{0} {1}/pathway-pathway.py -model {2} -net {3} -edge_node {3}\n".format(wfmo.get_python(),script_dir,model_network_dir,unetwork_dir))
                        os.system("{0} {1}/pathway-pathway.py -model {2} -net {3} -edge_node {3}".format(wfmo.get_python(), script_dir, model_network_dir, unetwork_dir))
                        # Down
                        olog.write("{0} {1}/pathway-pathway.py -model {2} -net {3} -edge_node {3}\n".format(wfmo.get_python(),script_dir,model_network_dir,dnetwork_dir))
                        os.system("{0} {1}/pathway-pathway.py -model {2} -net {3} -edge_node {3}".format(wfmo.get_python(), script_dir, model_network_dir, dnetwork_dir))

                    #DEP
                    if count_enrich_path("{}/KEGG/KEGG.enriched.txt".format(of)) > 0:
                        olog.write("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_all.txt {network_dir} all\n".format(
                                script_dir=script_dir, network_dir=network_dir))
                        os.system("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_all.txt {network_dir} all".format(
                                script_dir=script_dir, network_dir=network_dir))
                        olog.write("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all\n".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=network_dir))
                        os.system("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir,network_dir=network_dir))

                        olog.write("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_top5.txt {network_dir} part\n".format(
                                script_dir=script_dir, network_dir=network_dir))
                        os.system("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_top5.txt {network_dir} part".format(
                                script_dir=script_dir, network_dir=network_dir))

                        olog.write("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part\n".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=network_dir))
                        os.system("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir,network_dir=network_dir))
                    #UP
                    if count_enrich_path("{}/KEGG/KEGG.enriched.txt".format(ofu)) > 0:
                        olog.write("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_all.txt {network_dir} all\n".format(
                                script_dir=script_dir, network_dir=unetwork_dir))
                        os.system("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_all.txt {network_dir} all".format(
                                script_dir=script_dir, network_dir=unetwork_dir))
                        olog.write("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all\n".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=unetwork_dir))
                        os.system("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir,network_dir=unetwork_dir))

                        olog.write("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_top5.txt {network_dir} part\n".format(
                                script_dir=script_dir, network_dir=unetwork_dir))
                        os.system("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_top5.txt {network_dir} part".format(
                                script_dir=script_dir, network_dir=unetwork_dir))

                        olog.write("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part\n".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=unetwork_dir))
                        os.system("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir,network_dir=unetwork_dir))
                    #Down
                    if count_enrich_path("{}/KEGG/KEGG.enriched.txt".format(ofd)) > 0:
                        olog.write("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_all.txt {network_dir} all\n".format(
                                script_dir=script_dir, network_dir=dnetwork_dir))
                        os.system("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_all.txt {network_dir} all".format(
                                script_dir=script_dir, network_dir=dnetwork_dir))
                        olog.write("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all\n".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=dnetwork_dir))
                        os.system("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=dnetwork_dir))

                        olog.write("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_top5.txt {network_dir} part\n".format(
                                script_dir=script_dir, network_dir=dnetwork_dir))
                        os.system("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_top5.txt {network_dir} part".format(
                                script_dir=script_dir, network_dir=dnetwork_dir))

                        olog.write("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part\n".format(
                                python=wfmo.get_python(), script_dir=script_dir,model_network_dir=model_network_dir, network_dir=dnetwork_dir))
                        os.system("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=dnetwork_dir))

                    # 转换输入文件格式
                    olog.write("Rscript {0}/extractInput4network.R {1}/Input/Transformed_input.txt {1}/Input/Network_input.txt\n".format(bin, of, of))
                    os.system("Rscript {0}/extractInput4network.R {1}/Input/Transformed_input.txt {1}/Input/Network_input.txt".format(bin, of, of))
                    olog.write("Rscript {0}/extractInput4network.R {1}/Input/Transformed_input.txt {1}/Input/Network_input.txt\n".format(bin, ofu, ofu))
                    os.system("Rscript {0}/extractInput4network.R {1}/Input/Transformed_input.txt {1}/Input/Network_input.txt".format(bin, ofu, ofu))
                    olog.write("Rscript {0}/extractInput4network.R {1}/Input/Transformed_input.txt {1}/Input/Network_input.txt\n".format(bin, ofd, ofd))
                    os.system("Rscript {0}/extractInput4network.R {1}/Input/Transformed_input.txt {1}/Input/Network_input.txt".format(bin, ofd, ofd))

                    # DEP
                    olog.write("Rscript {script_dir}/pathway-gene-ppi.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} all\n".format(
                            script_dir=script_dir, network_dir=network_dir, ifold=of))
                    os.system("Rscript {script_dir}/pathway-gene-ppi.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} all".format(
                            script_dir=script_dir, network_dir=network_dir, ifold=of))

                    olog.write("Rscript {script_dir}/pathway-gene-ppi.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} all\n".format(
                            script_dir=script_dir, network_dir=unetwork_dir, ifold=ofu))
                    os.system("Rscript {script_dir}/pathway-gene-ppi.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} all".format(
                            script_dir=script_dir, network_dir=unetwork_dir, ifold=ofu))

                    olog.write("Rscript {script_dir}/pathway-gene-ppi.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} all\n".format(
                            script_dir=script_dir, network_dir=dnetwork_dir, ifold=ofd))
                    os.system("Rscript {script_dir}/pathway-gene-ppi.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} all".format(
                            script_dir=script_dir, network_dir=dnetwork_dir, ifold=ofd))

                    #DEP
                    olog.write("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all\n".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=network_dir))
                    os.system("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=network_dir))
                    #UP
                    olog.write("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all\n".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=unetwork_dir))
                    os.system("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=unetwork_dir))
                    #Down
                    olog.write("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all\n".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=dnetwork_dir))
                    os.system("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=dnetwork_dir))

                    #DEP
                    olog.write(
                        "Rscript {script_dir}/pathway-gene-ppi.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} part\n".format(
                            python=wfmo.get_python(), script_dir=script_dir, network_dir=network_dir, ifold=of))
                    os.system(
                        "Rscript {script_dir}/pathway-gene-ppi.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} part".format(
                            python=wfmo.get_python(), script_dir=script_dir, network_dir=network_dir, ifold=of))
                    #UP
                    olog.write(
                        "Rscript {script_dir}/pathway-gene-ppi.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} part\n".format(
                            python=wfmo.get_python(), script_dir=script_dir, network_dir=network_dir, ifold=of))
                    os.system(
                        "Rscript {script_dir}/pathway-gene-ppi.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} part".format(
                            python=wfmo.get_python(), script_dir=script_dir, network_dir=network_dir, ifold=of))
                    #Down
                    olog.write(
                        "Rscript {script_dir}/pathway-gene-ppi.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} part\n".format(
                            python=wfmo.get_python(), script_dir=script_dir, network_dir=network_dir, ifold=of))
                    os.system(
                        "Rscript {script_dir}/pathway-gene-ppi.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} part".format(
                            python=wfmo.get_python(), script_dir=script_dir, network_dir=network_dir, ifold=of))

                    #DEP
                    olog.write("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part\n".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=network_dir))
                    os.system("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=network_dir))
                    #UP
                    olog.write("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part\n".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=unetwork_dir))
                    os.system("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=unetwork_dir))
                    #Down
                    olog.write("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part\n".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=dnetwork_dir))
                    os.system("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=dnetwork_dir))

                    #DEP
                    os.system("Rscript {script_dir}/ppi.R {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} {symbol_idx} {fc_idx}".format(
                            script_dir=script_dir, ifold=of, network_dir=network_dir, symbol_idx=1, fc_idx=2))
                    #UP
                    os.system("Rscript {script_dir}/ppi.R {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} {symbol_idx} {fc_idx}".format(
                            script_dir=script_dir, ifold=ofu, network_dir=unetwork_dir, symbol_idx=1, fc_idx=2))
                    #Down
                    os.system("Rscript {script_dir}/ppi.R {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} {symbol_idx} {fc_idx}".format(
                            script_dir=script_dir, ifold=ofd, network_dir=dnetwork_dir, symbol_idx=1, fc_idx=2))


                    if os.path.exists(network_dir + "/node_edge/ppi_edge.txt") and os.path.exists(unetwork_dir + "/node_edge/ppi_edge.txt") and os.path.exists(
                            dnetwork_dir + "/node_edge/ppi_edge.txt"):
                        os.system('%s %s/ppi.py -model %s -net %s -edge_node %s' % (wfmo.get_python(), script_dir, model_network_dir, network_dir, network_dir))
                        os.system('%s %s/ppi.py -model %s -net %s -edge_node %s' % (wfmo.get_python(), script_dir, model_network_dir, unetwork_dir, unetwork_dir))
                        os.system('%s %s/ppi.py -model %s -net %s -edge_node %s' % (wfmo.get_python(), script_dir, model_network_dir, dnetwork_dir, dnetwork_dir))
                    print("\nNetwork plot done!\n")

                # 无fold change
                elif os.path.exists("{}/KEGG/KEGG.enriched.txt".format(of)) and os.path.exists("{}/KEGG/KEGG.enriched.txt".format(ofu)) and os.path.exists(
                        "{}/KEGG/KEGG.enriched.txt".format(ofd)):
                    olog.write("{0} {1}/get_path2gene.py -in {2}/KEGG/KEGG.enriched.txt -od {3}\n".format(wfmo.get_python(), bin, of, network_dir))
                    os.system("{0} {1}/get_path2gene.py -in {2}/KEGG/KEGG.enriched.txt -od {3}".format(wfmo.get_python(), bin, of, network_dir))

                    olog.write("{0} {1}/get_path2gene.py -in {2}/KEGG/KEGG.enriched.txt -od {3}\n".format(wfmo.get_python(), bin, ofu, unetwork_dir))
                    os.system("{0} {1}/get_path2gene.py -in {2}/KEGG/KEGG.enriched.txt -od {3}".format(wfmo.get_python(), bin, ofu, unetwork_dir))

                    olog.write("{0} {1}/get_path2gene.py -in {2}/KEGG/KEGG.enriched.txt -od {3}\n".format(wfmo.get_python(), bin, ofd, dnetwork_dir))
                    os.system("{0} {1}/get_path2gene.py -in {2}/KEGG/KEGG.enriched.txt -od {3}".format(wfmo.get_python(), bin, ofd, dnetwork_dir))

                    #DEP
                    if count_enrich_path("{}/KEGG/KEGG.enriched.txt".format(of)) > 0:
                        olog.write("Rscript {0}/pathway-pathway.R {1}/node_edge/edge_all.txt {1}\n".format(script_dir, network_dir))
                        os.system("Rscript {0}/pathway-pathway.R {1}/node_edge/edge_all.txt {1}".format(script_dir, network_dir))

                        olog.write("{0} {1}/pathway-pathway.py -model {2} -net {3} -edge_node {3}\n".format(wfmo.get_python(), script_dir, model_network_dir, network_dir))
                        os.system("{0} {1}/pathway-pathway.py -model {2} -net {3} -edge_node {3}".format(wfmo.get_python(), script_dir, model_network_dir, network_dir))

                        olog.write("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_all.txt {network_dir} all\n".format(
                                script_dir=script_dir, network_dir=network_dir))
                        os.system("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_all.txt {network_dir} all".format(
                                script_dir=script_dir, network_dir=network_dir))

                        olog.write("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all\n".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=network_dir))
                        os.system("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=network_dir))
                    #UP
                    if count_enrich_path("{}/KEGG/KEGG.enriched.txt".format(ofu)) > 0:
                        olog.write("Rscript {0}/pathway-pathway.R {1}/node_edge/edge_all.txt {1}\n".format(script_dir, unetwork_dir))
                        os.system("Rscript {0}/pathway-pathway.R {1}/node_edge/edge_all.txt {1}".format(script_dir, unetwork_dir))

                        olog.write("{0} {1}/pathway-pathway.py -model {2} -net {3} -edge_node {3}\n".format(wfmo.get_python(), script_dir, model_network_dir, unetwork_dir))
                        os.system("{0} {1}/pathway-pathway.py -model {2} -net {3} -edge_node {3}".format(wfmo.get_python(), script_dir, model_network_dir, unetwork_dir))

                        olog.write("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_all.txt {network_dir} all\n".format(
                                script_dir=script_dir, network_dir=unetwork_dir))
                        os.system("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_all.txt {network_dir} all".format(
                                script_dir=script_dir, network_dir=unetwork_dir))

                        olog.write("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all\n".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=unetwork_dir))
                        os.system("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=unetwork_dir))
                    #Down
                    if count_enrich_path("{}/KEGG/KEGG.enriched.txt".format(ofd)) > 0:
                        olog.write("Rscript {0}/pathway-pathway.R {1}/node_edge/edge_all.txt {1}\n".format(script_dir, dnetwork_dir))
                        os.system("Rscript {0}/pathway-pathway.R {1}/node_edge/edge_all.txt {1}".format(script_dir, dnetwork_dir))

                        olog.write("{0} {1}/pathway-pathway.py -model {2} -net {3} -edge_node {3}\n".format(wfmo.get_python(), script_dir, model_network_dir, dnetwork_dir))
                        os.system("{0} {1}/pathway-pathway.py -model {2} -net {3} -edge_node {3}".format(wfmo.get_python(), script_dir, model_network_dir, dnetwork_dir))

                        olog.write("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_all.txt {network_dir} all\n".format(
                                script_dir=script_dir, network_dir=dnetwork_dir))
                        os.system("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_all.txt {network_dir} all".format(
                                script_dir=script_dir, network_dir=dnetwork_dir))

                        olog.write("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all\n".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=dnetwork_dir))
                        os.system("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=dnetwork_dir))

                    #DEP
                    if count_enrich_path("{}/KEGG/KEGG.enriched.txt".format(of)) > 0:
                        olog.write("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_top5.txt {network_dir} part\n".format(
                                script_dir=script_dir, network_dir=network_dir))
                        os.system("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_top5.txt {network_dir} part".format(
                                script_dir=script_dir, network_dir=network_dir))

                        olog.write("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part\n".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=network_dir))
                        os.system("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=network_dir))
                    #UP
                    if count_enrich_path("{}/KEGG/KEGG.enriched.txt".format(ofu)) > 0:
                        olog.write("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_top5.txt {network_dir} part\n".format(
                                script_dir=script_dir, network_dir=unetwork_dir))
                        os.system("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_top5.txt {network_dir} part".format(
                                script_dir=script_dir, network_dir=unetwork_dir))

                        olog.write("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part\n".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=unetwork_dir))
                        os.system("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=unetwork_dir))
                    #Down
                    if count_enrich_path("{}/KEGG/KEGG.enriched.txt".format(ofd)) > 0:
                        olog.write("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_top5.txt {network_dir} part\n".format(
                                script_dir=script_dir, network_dir=dnetwork_dir))
                        os.system("Rscript {script_dir}/pathway-gene.R {network_dir}/node_edge/edge_top5.txt {network_dir} part".format(
                                script_dir=script_dir, network_dir=dnetwork_dir))

                        olog.write("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part\n".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=dnetwork_dir))
                        os.system("{python} {script_dir}/pathway-gene.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part".format(
                                python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=dnetwork_dir))

                    #DEP
                    olog.write("Rscript {script_dir}/pathway-gene-ppi_fc1.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {network_dir} all\n".format(
                            script_dir=script_dir, network_dir=network_dir, ifold=of))
                    os.system("Rscript {script_dir}/pathway-gene-ppi_fc1.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {network_dir} all".format(
                            script_dir=script_dir, network_dir=network_dir, ifold=of))
                    #UP
                    olog.write("Rscript {script_dir}/pathway-gene-ppi_fc1.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {network_dir} all\n".format(
                            script_dir=script_dir, network_dir=unetwork_dir, ifold=ofu))
                    os.system("Rscript {script_dir}/pathway-gene-ppi_fc1.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {network_dir} all".format(
                            script_dir=script_dir, network_dir=unetwork_dir, ifold=ofu))
                    #Down
                    olog.write("Rscript {script_dir}/pathway-gene-ppi_fc1.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {network_dir} all\n".format(
                            script_dir=script_dir, network_dir=dnetwork_dir, ifold=ofd))
                    os.system("Rscript {script_dir}/pathway-gene-ppi_fc1.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {network_dir} all".format(
                            script_dir=script_dir, network_dir=dnetwork_dir, ifold=ofd))

                    #DEP
                    olog.write("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all\n".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=network_dir))
                    os.system("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=network_dir))
                    #UP
                    olog.write("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all\n".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=unetwork_dir))
                    os.system("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=unetwork_dir))
                    #Down
                    olog.write("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all\n".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=dnetwork_dir))
                    os.system("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s all".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=dnetwork_dir))

                    #DEP
                    olog.write(
                        "Rscript {script_dir}/pathway-gene-ppi.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {ifold}/DESelection/DEPProtein.xls {network_dir} part\n".format(
                            python=wfmo.get_python(), script_dir=script_dir, network_dir=network_dir, ifold=of))
                    os.system(
                        "Rscript {script_dir}/pathway-gene-ppi.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {ifold}/DESelection/DEPProtein.xls {network_dir} part".format(
                            python=wfmo.get_python(), script_dir=script_dir, network_dir=network_dir, ifold=of))
                    #UP
                    olog.write(
                        "Rscript {script_dir}/pathway-gene-ppi.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {ifold}/DESelection/DEPProtein.xls {network_dir} part\n".format(
                            python=wfmo.get_python(), script_dir=script_dir, network_dir=unetwork_dir, ifold=ofu))
                    os.system(
                        "Rscript {script_dir}/pathway-gene-ppi.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {ifold}/DESelection/DEPProtein.xls {network_dir} part".format(
                            python=wfmo.get_python(), script_dir=script_dir, network_dir=unetwork_dir, ifold=ofu))
                    #Down
                    olog.write(
                        "Rscript {script_dir}/pathway-gene-ppi.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {ifold}/DESelection/DEPProtein.xls {network_dir} part\n".format(
                            python=wfmo.get_python(), script_dir=script_dir, network_dir=dnetwork_dir, ifold=ofd))
                    os.system(
                        "Rscript {script_dir}/pathway-gene-ppi.R {network_dir}/node_edge/edge_top5.txt {ifold}/PPI/NonRPPI.txt {ifold}/DESelection/DEPProtein.xls {network_dir} part".format(
                            python=wfmo.get_python(), script_dir=script_dir, network_dir=dnetwork_dir, ifold=ofd))

                    #DEP
                    olog.write("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part\n".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=network_dir))
                    os.system("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=network_dir))
                    #UP
                    olog.write("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part\n".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=unetwork_dir))
                    os.system("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=unetwork_dir))
                    #Down
                    olog.write("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part\n".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=dnetwork_dir))
                    os.system("{python} {script_dir}/pathway-gene-ppi.py -model {model_network_dir} -net {network_dir} -edge_node {network_dir} -s part".format(
                            python=wfmo.get_python(), script_dir=script_dir, model_network_dir=model_network_dir, network_dir=dnetwork_dir))

                    # 需要构造一个FC信息的文件
                    construct_fc_file("{}/Input/Transformed_input.txt".format(of), "{}/Input/Network_input.txt".format(of))
                    construct_fc_file("{}/Input/Transformed_input.txt".format(ofu), "{}/Input/Network_input.txt".format(ofu))
                    construct_fc_file("{}/Input/Transformed_input.txt".format(ofd), "{}/Input/Network_input.txt".format(ofd))

                    #DEP
                    olog.write("Rscript {script_dir}/ppi.R {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} {symbol_idx} {fc_idx}\n".format(
                            script_dir=script_dir, ifold=of, network_dir=network_dir, symbol_idx=1, fc_idx=2))
                    os.system("Rscript {script_dir}/ppi.R {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} {symbol_idx} {fc_idx}".format(
                            script_dir=script_dir, ifold=of, network_dir=network_dir, symbol_idx=1, fc_idx=2))
                    #UP
                    olog.write("Rscript {script_dir}/ppi.R {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} {symbol_idx} {fc_idx}\n".format(
                            script_dir=script_dir, ifold=ofu, network_dir=unetwork_dir, symbol_idx=1, fc_idx=2))
                    os.system("Rscript {script_dir}/ppi.R {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} {symbol_idx} {fc_idx}".format(
                            script_dir=script_dir, ifold=ofu, network_dir=unetwork_dir, symbol_idx=1, fc_idx=2))
                    #Down
                    olog.write("Rscript {script_dir}/ppi.R {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} {symbol_idx} {fc_idx}\n".format(
                            script_dir=script_dir, ifold=ofd, network_dir=dnetwork_dir, symbol_idx=1, fc_idx=2))
                    os.system("Rscript {script_dir}/ppi.R {ifold}/PPI/NonRPPI.txt {ifold}/Input/Network_input.txt {network_dir} {symbol_idx} {fc_idx}".format(
                            script_dir=script_dir, ifold=ofd, network_dir=dnetwork_dir, symbol_idx=1, fc_idx=2))

                    if os.path.exists(network_dir + "/node_edge/ppi_edge.txt") and os.path.exists(unetwork_dir + "/node_edge/ppi_edge.txt") and os.path.exists(
                            dnetwork_dir + "/node_edge/ppi_edge.txt"):
                        olog.write('%s %s/ppi.py -model %s -net %s -edge_node %s\n' % (wfmo.get_python(), script_dir, model_network_dir, network_dir, network_dir))
                        os.system('%s %s/ppi.py -model %s -net %s -edge_node %s' % (wfmo.get_python(), script_dir, model_network_dir, network_dir, network_dir))
                        #UP
                        olog.write('%s %s/ppi.py -model %s -net %s -edge_node %s\n' % (wfmo.get_python(), script_dir, model_network_dir, unetwork_dir, unetwork_dir))
                        os.system('%s %s/ppi.py -model %s -net %s -edge_node %s' % (wfmo.get_python(), script_dir, model_network_dir, unetwork_dir, unetwork_dir))
                        #Down
                        olog.write('%s %s/ppi.py -model %s -net %s -edge_node %s\n' % (wfmo.get_python(), script_dir, model_network_dir, dnetwork_dir, dnetwork_dir))
                        os.system('%s %s/ppi.py -model %s -net %s -edge_node %s' % (wfmo.get_python(), script_dir, model_network_dir, dnetwork_dir, dnetwork_dir))



"""多组比较数据整合绘制柱状图"""
def merge_updown_bar(L, pn, bin, olog):
    if len(L) > 0:
        pass
        s = ''
        for comp in L:
            s = s + ' ' + pn + '/' + comp + '/DESelection/updown.stat.txt ' + comp
        # print("Rscript {}/mergeUpdownStat.R {} {}".format(bin, pn, s))
        olog.write("Rscript {}/mergeUpdownStat.R {} {}\n".format(bin, pn, s))
        os.system("Rscript {}/mergeUpdownStat.R {} {}".format(bin, pn, s))

        if os.path.exists(pn + '/multiUpdown.stat.txt'):
            olog.write("Rscript {}/dUpdownbar.R {}/multiUpdown.stat.txt {}\n".format(bin, pn, pn))
            os.system("Rscript {}/dUpdownbar.R {}/multiUpdown.stat.txt {}".format(bin, pn, pn))
"""维恩图"""
def venn_anal(L, pn, bin, olog):

    if len(L) > 1 and len(L) <= 5:
        L = sorted(L)
        s = ''
        sup = ''
        sdown = ''
        for comp in L:
            s = s + ' ' + pn + '/' + comp + '/DESelection/DEP.xls'
            sup = sup + ' ' + pn + '/' + comp + '/DESelection/Up.xls'
            sdown = sdown + ' ' + pn + '/' + comp + '/DESelection/Down.xls'

        if not os.path.exists(pn + '/Venn'):
            os.makedirs(pn + '/Venn')

        olog.write("Rscript {}/dVenn.R {} {} \"{}\"\n".format(bin, s, pn + '/Venn/all', ';'.join(L)))
        os.system("Rscript {}/dVenn.R {} {} \"{}\"".format(bin, s, pn + '/Venn/all', ';'.join(L)))

        olog.write("Rscript {}/dVenn.R {} {} \"{}\"\n".format(bin, sup, pn + '/Venn/up', ';'.join(L)))
        os.system("Rscript {}/dVenn.R {} {} \"{}\"".format(bin, sup, pn + '/Venn/up', ';'.join(L)))

        olog.write("Rscript {}/dVenn.R {} {} \"{}\"\n".format(bin, sdown, pn + '/Venn/down', ';'.join(L)))
        os.system("Rscript {}/dVenn.R {} {} \"{}\"".format(bin, sdown, pn + '/Venn/down', ';'.join(L)))


def count_enrich_path(infile, pvalue_idx=4, pvalue_cutoff=0.05):
    ifile = open(infile, 'r')
    c = 0
    for line in ifile.readlines()[1:]:
        row = line.rstrip().split('\t')
        if float(row[pvalue_idx]) < float(pvalue_cutoff):
            c = c + 1
    ifile.close()
    return c

def construct_fc_file(infile, outfile):
    if os.path.exists(infile):
        ih = open(infile)
        oh = open(outfile,'w')
        lines = ih.readlines()
        title=lines[0]
        titlelist = title.rstrip().split('\t')
        oh.write("{}\t{}\n".format(titlelist[1], 'Log2FC'))
        for line in lines[1:]:
            row = line.rstrip().split('\t')
            oh.write("{}\t{}\n".format(row[1], 1))

        ih.close()
        oh.close()
    else:
        print("找不到输入文件")

"""挑选文件到结果目录"""
def file_selection(peptideFile, L, pn, olog):
    reportFold = pn + '/' + pn + '_Result'
    if os.path.exists(reportFold):
        print("结果目录已存在")
        sys.exit()

    if not os.path.exists(reportFold):
        os.makedirs(reportFold)

    createdir(reportFold + '/1.ProjectInfo')

    # result1
    renamefile(peptideFile, reportFold + '/1.ProjectInfo/Peptide.xlsx')
    cpfile(pn + '/1.Info/ProbabilitySite.xlsx', reportFold + '/1.ProjectInfo')
    cpfile(pn + '/1.Info/DEPProbabilitySite.xlsx', reportFold + '/1.ProjectInfo')
    cpfile(pn + '/1.Info/Protein_PTM.xlsx', reportFold + '/1.ProjectInfo')

    # result
    """
    if os.path.exists(pn + '/2.QualityControl') and not os.path.exists(reportFold + '/2.QualityControl'):
        shutil.copytree(pn + '/2.QualityControl', reportFold + '/2.QualityControl')
        if os.path.exists(reportFold + '/2.QualityControl/peptide_error.png'):
            os.remove(reportFold + '/2.QualityControl/peptide_error.png')
        if os.path.exists(reportFold + '/2.QualityControl/peptide_error.pdf'):
            os.remove(reportFold + '/2.QualityControl/peptide_error.pdf')
    """
    # result2
    if os.path.exists(pn + '/2.SampleAnalysis') and not os.path.exists(reportFold + '/2.SampleAnalysis'):
        shutil.copytree(pn + '/2.SampleAnalysis', reportFold + '/2.SampleAnalysis')


    # result3 总蛋白结果整理
    allProFold = reportFold + '/3.AllProtein'
    createdir(allProFold)
    cpfile(pn + '/summary.xlsx', allProFold)

    createdir(allProFold + '/3.1GO')
    goEnrichFold = allProFold + '/3.1GO'
    # go level2 file
    cpfile(pn + '/AllProtein/GOLevel2/GOLevel2.xlsx', goEnrichFold)
    cpfile(pn + '/AllProtein/GOLevel2/GOlevel2CountBar.png', goEnrichFold)
    cpfile(pn + '/AllProtein/GOLevel2/GOlevel2CountBar.pdf', goEnrichFold)

    # go enrich file

    if os.path.exists(pn + '/AllProtein/GOEnrich_Species'):
        cpfile(pn + '/AllProtein/GOEnrich_Species/pCountPoint.pdf', goEnrichFold)
        cpfile(pn + '/AllProtein/GOEnrich_Species/pCountPoint.png', goEnrichFold)
        cpfile(pn + '/AllProtein/GOEnrich_Species/pRFPoint.pdf', goEnrichFold)
        cpfile(pn + '/AllProtein/GOEnrich_Species/pRFPoint.png', goEnrichFold)

        cpfile(pn + '/AllProtein/GOEnrich_Species/BP.EnrichedBar.pdf', goEnrichFold)
        cpfile(pn + '/AllProtein/GOEnrich_Species/BP.EnrichedBar.png', goEnrichFold)
        cpfile(pn + '/AllProtein/GOEnrich_Species/MF.EnrichedBar.pdf', goEnrichFold)
        cpfile(pn + '/AllProtein/GOEnrich_Species/MF.EnrichedBar.png', goEnrichFold)
        cpfile(pn + '/AllProtein/GOEnrich_Species/CC.EnrichedBar.pdf', goEnrichFold)
        cpfile(pn + '/AllProtein/GOEnrich_Species/CC.EnrichedBar.png', goEnrichFold)

        cpfile(pn + '/AllProtein/GOEnrich_Species/BP.DAG.pdf', goEnrichFold)
        cpfile(pn + '/AllProtein/GOEnrich_Species/BP.DAG.png', goEnrichFold)
        cpfile(pn + '/AllProtein/GOEnrich_Species/MF.DAG.pdf', goEnrichFold)
        cpfile(pn + '/AllProtein/GOEnrich_Species/MF.DAG.png', goEnrichFold)
        cpfile(pn + '/AllProtein/GOEnrich_Species/CC.DAG.pdf', goEnrichFold)
        cpfile(pn + '/AllProtein/GOEnrich_Species/CC.DAG.png', goEnrichFold)

        renamefile(pn + '/AllProtein/GOEnrich_Species/Top10BP.EnrichedSymbol.pdf',
                   goEnrichFold + '/BP.RichFactor.pdf')
        renamefile(pn + '/AllProtein/GOEnrich_Species/Top10BP.EnrichedSymbol.png',
                   goEnrichFold + '/BP.RichFactor.png')
        renamefile(pn + '/AllProtein/GOEnrich_Species/Top10MF.EnrichedSymbol.pdf',
                   goEnrichFold + '/MF.RichFactor.pdf')
        renamefile(pn + '/AllProtein/GOEnrich_Species/Top10MF.EnrichedSymbol.png',
                   goEnrichFold + '/MF.RichFactor.png')
        renamefile(pn + '/AllProtein/GOEnrich_Species/Top10CC.EnrichedSymbol.pdf',
                   goEnrichFold + '/CC.RichFactor.pdf')
        renamefile(pn + '/AllProtein/GOEnrich_Species/Top10CC.EnrichedSymbol.png',
                   goEnrichFold + '/CC.RichFactor.png')

        renamefile(pn + '/AllProtein/GOEnrich_Species/Top20.BP.EnrichedSymbol2.pdf',
                   goEnrichFold + '/BP.EnrichedSymbol.pdf')
        renamefile(pn + '/AllProtein/GOEnrich_Species/Top20.BP.EnrichedSymbol2.png',
                   goEnrichFold + '/BP.EnrichedSymbol.png')
        renamefile(pn + '/AllProtein/GOEnrich_Species/Top20.MF.EnrichedSymbol2.pdf',
                   goEnrichFold + '/MF.EnrichedSymbol.pdf')
        renamefile(pn + '/AllProtein/GOEnrich_Species/Top20.MF.EnrichedSymbol2.png',
                   goEnrichFold + '/MF.EnrichedSymbol.png')
        renamefile(pn + '/AllProtein/GOEnrich_Species/Top20.CC.EnrichedSymbol2.pdf',
                   goEnrichFold + '/CC.EnrichedSymbol.pdf')
        renamefile(pn + '/AllProtein/GOEnrich_Species/Top20.CC.EnrichedSymbol2.png',
                   goEnrichFold + '/CC.EnrichedSymbol.png')

        cpfile(pn + '/AllProtein/GOEnrich_Species/GOEnrich.xlsx', goEnrichFold)

    # KEGG files
    allProtKegg = pn + '/AllProtein/KEGG'

    keggfold = allProFold + '/3.2KEGG'
    createdir(keggfold)
    if os.path.exists(allProtKegg):
        cpfile(allProtKegg + '/KEGG.Sig.Bar.png', keggfold)
        cpfile(allProtKegg + '/KEGG.Sig.Bar.pdf', keggfold)

        cpfile(allProtKegg + '/KEGG.Sig.Bar2.png', keggfold)
        cpfile(allProtKegg + '/KEGG.Sig.Bar2.pdf', keggfold)

        cpfile(allProtKegg + '/Top10EnrichedBar.png', keggfold)
        cpfile(allProtKegg + '/Top10EnrichedBar.pdf', keggfold)

        cpfile(allProtKegg + '/Top10EnrichSymbol.png', keggfold)
        cpfile(allProtKegg + '/Top10EnrichSymbol.pdf', keggfold)

        cpfile(allProtKegg + '/Top20EnrichedSymbol2.png', keggfold)
        cpfile(allProtKegg + '/Top20EnrichedSymbol2.pdf', keggfold)

        renamefile(allProtKegg + '/KEGG.bubble.png', keggfold + '/kegg.sig.png')
        renamefile(allProtKegg + '/KEGG.bubble.pdf', keggfold + '/kegg.sig.pdf')

        cpfile(allProtKegg + '/KEGG.enriched.xlsx', keggfold)

    sublocfold = allProFold + '/3.3subcellular_localization'
    # createdir(sublocfold)
    if os.path.exists(pn + '/AllProtein/GO_subcellular_localization'):
        shutil.copytree(pn + '/AllProtein/GO_subcellular_localization', sublocfold)

    ##### result4 所有可信位点Motif分析
    allMotifFold = reportFold + '/4.MotifAnalysis'
    shutil.copytree(pn + '/MotifAnalysis', allMotifFold)

    ####result5 差异筛选结果整理
    depstatfold = reportFold + '/5.DEPStatistics'
    createdir(depstatfold)

    cpfile(pn + '/multiDEPStatbar.png', depstatfold)
    cpfile(pn + '/multiDEPStatbar.pdf', depstatfold)

    # venn diagram
    vennfold = depstatfold + '/Venn'
    if os.path.exists(pn + '/Venn/all.tif'):
        createdir(vennfold)
        venntype = ['all', 'up', 'down']
        for vt in venntype:
            cpfile(pn + '/Venn/{}.tif'.format(vt), vennfold)
            cpfile(pn + '/Venn/{}.merge.xls'.format(vt), vennfold)
            # cpfile(pn + '/Venn/{}.euler.png'.format(vt), vennfold)
            # cpfile(pn + '/Venn/{}.euler.pdf'.format(vt), vennfold)
            # cpfile(pn + '/Venn/{}.eulerEllipse.png'.format(vt), vennfold)
            # cpfile(pn + '/Venn/{}.eulerEllipse.pdf'.format(vt), vennfold)
            # cpfile(pn + '/Venn/{}.CommonAcc.xls'.format(vt), vennfold)

    ####result7 差异功能分析结果
    depfuncfold = reportFold + '/6.DEPFunction'
    createdir(depfuncfold)

    # 拷贝图片
    for comp in L:
        createdir(depfuncfold + '/' + comp)
        createdir(depstatfold + '/' + comp)
        DESelection = pn + '/' + comp + '/DESelection'
        renamefile(DESelection + '/Volcano.pdf', depstatfold + '/' + comp + '/' + comp + '.volcano.pdf')
        renamefile(DESelection + '/Volcano.png', depstatfold + '/' + comp + '/' + comp + '.volcano.png')

        renamefile(DESelection + '/Volcano_style2.pdf', depstatfold + '/' + comp + '/' + comp + '.volcano_s2.pdf')
        renamefile(DESelection + '/Volcano_style2.png', depstatfold + '/' + comp + '/' + comp + '.volcano_s2.png')

        renamefile(DESelection + '/FoldChange.scatter.pdf',
                   depstatfold + '/' + comp + '/' + comp + '.FoldChange.scatter.pdf')
        renamefile(DESelection + '/FoldChange.scatter.png',
                   depstatfold + '/' + comp + '/' + comp + '.FoldChange.scatter.png')

        renamefile(DESelection + '/pvalue.scatter.pdf', depstatfold + '/' + comp + '/' + comp + '.pvalue.scatter.pdf')
        renamefile(DESelection + '/pvalue.scatter.png', depstatfold + '/' + comp + '/' + comp + '.pvalue.scatter.png')

        renamefile(DESelection + '/updown.bar.pdf', depstatfold + '/' + comp + '/' + comp + '.updown.bar.pdf')
        renamefile(DESelection + '/updown.bar.png', depstatfold + '/' + comp + '/' + comp + '.updown.bar.png')

        renamefile(DESelection + '/DEP.Donut.pdf', depstatfold + '/' + comp + '/' + comp + '.DEP.Donut.pdf')
        renamefile(DESelection + '/DEP.Donut.png', depstatfold + '/' + comp + '/' + comp + '.DEP.Donut.png')

        renamefile(DESelection + '/Accession.pheatmap.pdf',
                   depstatfold + '/' + comp + '/' + comp + '.accession.pheatmap.pdf')
        renamefile(DESelection + '/Accession.pheatmap.png',
                   depstatfold + '/' + comp + '/' + comp + '.accession.pheatmap.png')

        renamefile(DESelection + '/Accession.pheatmapcolcluster.pdf',
                   depstatfold + '/' + comp + '/' + comp + '.accession.pheatmapcolcluster.pdf')
        renamefile(DESelection + '/Accession.pheatmapcolcluster.png',
                   depstatfold + '/' + comp + '/' + comp + '.accession.pheatmapcolcluster.png')

        # renamefile(DESelection + '/genesymbol.pheatmap.pdf',
        #            depstatfold + '/' + comp + '/' + comp + '.genesymbol.pheatmap.pdf')
        # renamefile(DESelection + '/genesymbol.pheatmap.png',
        #            depstatfold + '/' + comp + '/' + comp + '.genesymbol.pheatmap.png')

        # renamefile(DESelection + '/genesymbol.pheatmap.pdf',
        #            depstatfold + '/' + comp + '/' + comp + '.genesymbol.pheatmapcolcluster.pdf')
        # renamefile(DESelection + '/genesymbol.pheatmap.png',
        #            depstatfold + '/' + comp + '/' + comp + '.genesymbol.pheatmapcolcluster.png')


        # DEP单个绘图的拷贝
        shutil.copytree(DESelection + '/Jitter', depstatfold + '/' + comp + '/' + comp + '.jitter')
        shutil.copytree(DESelection + '/Boxplot', depstatfold + '/' + comp + '/' + comp + '.box')
        shutil.copytree(DESelection + '/ScatterSE', depstatfold + '/' + comp + '/' + comp + '.scatter')
        shutil.copytree(DESelection + '/BarComp', depstatfold + '/' + comp + '/' + comp + '.bar')

        # 功能分析结果
        renamefile(pn + '/' + comp + '/all_summary.xlsx', depfuncfold + '/' + comp + '/' + comp + '.all_summary.xlsx')
        renamefile(pn + '/' + comp + '/up_summary.xlsx', depfuncfold + '/' + comp + '/' + comp + '.up_summary.xlsx')
        renamefile(pn + '/' + comp + '/down_summary.xlsx', depfuncfold + '/' + comp + '/' + comp + '.down_summary.xlsx')

        # 功能分析结果
        for type in ['all', 'up', 'down']:
            # # 以全部鉴定到的蛋白为参照的结果
            # goidentIn = pn + '/' + comp + '/' + type + '/GOEnrich_Identified'
            # if os.path.exists(goidentIn):
            #     goidentOut = depfuncfold + '/' + comp + '/' + comp + '.' + type + '_GO_Ident'
            #     shutil.copytree(goidentIn, goidentOut)
            #     deleteFiles = ['BP.addFC.txt', 'BP.cnet.pdf', 'BP.txt', 'BP.updown.txt',
            #                    'MF.addFC.txt', 'MF.cnet.pdf', 'MF.txt', 'MF.updown.txt',
            #                    'CC.addFC.txt', 'CC.cnet.pdf', 'CC.txt', 'CC.updown.txt']
            #     for file in deleteFiles:
            #         if os.path.exists(goidentOut + '/' + file):
            #             os.remove(goidentOut + '/' + file)

            # 以物种全蛋白为参照的结果
            gospeciesIn = pn + '/' + comp + '/' + type + '/GOEnrich_Species'
            gospeciesOut = depfuncfold + '/' + comp + '/' + comp + '.' + type + '_GO_Species'
            if os.path.exists(gospeciesIn):
                shutil.copytree(gospeciesIn, gospeciesOut)
                deleteFiles = ['BP.addFC.txt', 'BP.cnet.pdf', 'BP.txt', 'BP.updown.txt',
                               'MF.addFC.txt', 'MF.cnet.pdf', 'MF.txt', 'MF.updown.txt',
                               'CC.addFC.txt', 'CC.cnet.pdf', 'CC.txt', 'CC.updown.txt']
                for file in deleteFiles:
                    if os.path.exists(gospeciesOut + '/' + file):
                        os.remove(gospeciesOut + '/' + file)
            # GO level2
            golevelIn = pn + '/' + comp + '/' + type + '/GOLevel2'
            if os.path.exists(golevelIn):
                shutil.copytree(golevelIn, depfuncfold + '/' + comp + '/' + comp + '.' + type + '_GO_Level2')

            #亚细胞定位
            gocelllocIn = pn + '/' + comp + '/' + type + '/GO_subcellular_localization'
            gocelllocOut =depfuncfold + '/' + comp +'/'+comp+'.'+type+'_subcellular_localization'
            if os.path.exists(gocelllocIn):
                shutil.copytree(gocelllocIn,gocelllocOut)

            # KEGG
            keggIn = pn + '/' + comp + '/' + type + '/KEGG'
            keggOut = depfuncfold + '/' + comp + '/' + comp + '.' + type + '_KEGG'
            if os.path.exists(keggIn):
                shutil.copytree(keggIn, keggOut)
                deleteFiles = ['KEGG.enriched.xls', 'KEGG.enriched.txt']
                for file in deleteFiles:
                    if os.path.exists(keggOut + '/' + file):
                        os.remove(keggOut + '/' + file)

            # if type == 'all':
            #     """
            #     cpfile(keggOut + '/Circos/circos.svg', keggOut)
            #     cpfile(keggOut + '/Circos/circos.png', keggOut)
            #     cpfile(keggOut + '/Circos/legend.svg', keggOut)
            #     """
            #     if os.path.exists(keggOut + '/Circos'):
            #         shutil.rmtree(keggOut + '/Circos')

            pfin_in = pn + '/' + comp + '/' + type + '/Network'
            pfin_out = depfuncfold + '/' + comp + '/' + comp + '.' + type + '_PFIN'
            if os.path.exists(pfin_in):
                shutil.copytree(pfin_in, pfin_out)
            else:
                print(pfin_in + "路径不存在！\n")

        #拷贝KEGG updownSbar
        if os.path.exists(pn+'/'+comp+'/updownSbar.png'):
            cpfile(pn+'/'+comp+'/updownSbar.png', depfuncfold + '/' + comp + '/' + comp + '.' + 'all_KEGG')
            cpfile(pn+'/'+comp+'/updownSbar.pdf', depfuncfold + '/' + comp + '/' + comp + '.' + 'all_KEGG')

    ####result7 差异Motif分析结果
    depmotiffold = reportFold + '/7.DEPMotif'
    createdir(depmotiffold)
    for comp in L:
        for type in ['All', 'Up', 'Down']:
            # 差异Motif分析结果
            depmotifIn = pn + '/DEPMotif/' + comp + '/' + type
            depmotifOut = depmotiffold + '/' + comp + '/' + comp + '.' + type
            if os.path.exists(depmotifIn):
                shutil.copytree(depmotifIn, depmotifOut)
            else:
                print(depmotifIn + "路径不存在！\n")

#### Anova分析 ???
    if os.path.exists(pn + '/Anova'):
        shutil.copytree(pn + '/Anova', reportFold + '/7.ANOVA')
        delete_files = ['group.txt', 'input.txt',
                        'GOEnrich_Species/BP.txt', 'GOEnrich_Species/CC.txt','GOEnrich_Species/MF.txt',
                        'KEGG/KEGG.enriched.txt','KEGG/KEGG.enriched.xls']
        delete_folds = ['ReactomePA_Species', 'WikiPathway']

        ag = os.listdir(reportFold + '/7.ANOVA')
        for g in ag:
            delete_filelist(reportFold + '/7.ANOVA/'+ g, delete_files)
            delete_foldlist(reportFold + '/7.ANOVA/'+ g, delete_folds)

        for g in ag:
            filepath = reportFold + '/7.ANOVA/' + g + '/kmeans/'
            kmeansfile = filepath + 'zscore10Hartigan-Wong.Class.xls'
            if os.path.exists(kmeansfile):
                os.rename(kmeansfile, filepath + 'zscoreHartigan-Wong.Class.xls')

            kmeansfile = filepath + 'zscore10Hartigan-Wong.kmeans.png'
            if os.path.exists(kmeansfile):
                os.rename(kmeansfile, filepath + 'zscoreHartigan-Wong.kmeans.png')

            kmeansfile = filepath + 'zscore10Hartigan-Wong.kmeans.pdf'
            if os.path.exists(kmeansfile):
                os.rename(kmeansfile, filepath + 'zscoreHartigan-Wong.kmeans.pdf')

            kmeansfile = filepath + 'zscore10Hartigan-Wong.kmeansldf.csv'
            if os.path.exists(kmeansfile):
                os.rename(kmeansfile, filepath + 'zscoreHartigan-Wong.kmeansldf.csv')

            kmeansfile = filepath + 'zscoreavgk5Hartigan-Wong.Class.xls'
            if os.path.exists(kmeansfile):
                os.rename(kmeansfile, filepath + 'zscoreavgHartigan-Wong.Class.xls')

            kmeansfile = filepath + 'zscoreavgk5Hartigan-Wong.kmeans.png'
            if os.path.exists(kmeansfile):
                os.rename(kmeansfile, filepath + 'zscoreavgHartigan-Wong.kmeans.png')

            kmeansfile = filepath + 'zscoreavgk5Hartigan-Wong.kmeans.pdf'
            if os.path.exists(kmeansfile):
                os.rename(kmeansfile, filepath + 'zscoreavgHartigan-Wong.kmeans.pdf')

            kmeansfile = filepath + 'zscoreavgk5Hartigan-Wong.kmeansldf.csv'
            if os.path.exists(kmeansfile):
                os.rename(kmeansfile, filepath + 'zscoreavgHartigan-Wong.kmeansldf.csv')

"""按文件列表删除文件"""
def delete_filelist(fold, filelist):
    for f in filelist:
        if os.path.exists(fold + '/' + f):
            os.remove(fold + '/' + f)
"""按目录列表删除目录"""
def delete_foldlist(fold, foldlist):
    for f in foldlist:
        if os.path.exists(fold + '/' + f):
            shutil.rmtree(fold + '/' + f)

"""转换ID"""
def convert_identifier(infile, org, type, pn, bin, olog):
    olog.write("Rscript {}/convert_identifier.R {} {} {} {}\n".format(bin, infile, pn, org, type))
    os.system("Rscript {}/convert_identifier.R {} {} {} {}".format(bin, infile, pn, org, type))

"""创建一个目录"""
def createdir(d):
    if not os.path.exists(d):
        os.makedirs(d)

"""拷贝文件到目录"""
def cpfile(file, d):
    if os.path.exists(file) and os.path.exists(d):
        shutil.copy(file, d)
    else:
        print("Copy {} to {} failed".format(file, d))

"""拷贝文件到另一个文件名"""
def renamefile(file, tofile):
    if os.path.exists(file):
        shutil.copyfile(file, tofile)
    else:
        print("Cannot find {}".format(file))


"""获取当前程序所在的目录"""
def get_absdir():
    if not sys.argv[0].endswith('.py') and not sys.argv[0].endswith('.pyw'):
        absDir = os.path.dirname(sys.executable)
    else:
        absDir = os.path.dirname(os.path.abspath(__file__))

    return absDir


def get_bin():
    if platform.system() == 'Windows':
        bin = get_absdir()
    else:
        bin = '/usr/local/bio/PAA/bin'
    return bin


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', action='store', dest='pn', default='', help='Project ID_fc?? as output fold')
    parser.add_argument('-pf', action='store', dest='pf', default='Proteins.xlsx', help='input Protein.xlsx')
    parser.add_argument('-pef', action='store', dest='pef', default='PeptidesGroup.xlsx', help='input PeptideGroup.xlsx')
    parser.add_argument('-psm', action='store', dest='psm', default='PSMs.xlsx', help='input PSMs.xlsx')
    parser.add_argument('-info', action='store', dest='info', default='parameter.xlsx', help='input parameters.xlsx')
    parser.add_argument('-fc', action='store', dest='fc', default=1.5, help='fold change cutoff, default is 1.5')
    parser.add_argument('-pvalue', action='store', dest='pvalue', default=0.05, help='p-value cutoff, default is 0.05')
    parser.add_argument('-col', action='store', dest='col', default=3, help='columns number of facet drawing,default 3')
    parser.add_argument('-sm', action='store', dest='scaleMethod', default='log2', help='scale method for data draw')
    parser.add_argument('-org', action='store', dest='org', default='hsa', help='org code')
    parser.add_argument('-celllocation', action='store', dest='celllocation', default=1,
                        help='1:animal; 2:plant; 3:bacterial; 4:archaeal')
    parser.add_argument('-taxid', action='store', dest='taxid', default=9606, help='taxonomy ID, default is 9606 for Human')
    parser.add_argument('-reporttype', action='store', dest='reporttype', default='DIA', help="报告类型，TMT|LFQ|DIA, 默认是DIA")
    parser.add_argument('-reanal', action='store', dest='reanal', default='', help="重做: ppi|select|report")
    parser.add_argument('-mm', action='store', dest='maskminvalue', default=False,
                        help="是否屏蔽掉最小值后再做箱线图，默认是False, 通常Label要设为True, DIA看处理前是不是有缺失值")
    parser.add_argument('-isft', action='store', dest='isft', default='s', help='Search the library software. s:Spectronaut; m:MSFragger; p:PD')
    parser.add_argument('-ibg', action='store', dest='ibg', default='none', help='background.input file')
    parser.add_argument('-type', action='store', dest='type', default='symbol', help='convert type: symbol, accsymbol,accsymbolfc, accession')


    supp = ['hsa', 'rno', 'mmu', 'ssc', 'gga', 'bta', 'cel', 'ath', 'dre', 'dme', 'sce']
    code2taxid = dict()
    code2taxid['hsa'] = 9606
    code2taxid['rno'] = 10116
    code2taxid['mmu'] = 10090
    code2taxid['ssc'] = 9823
    code2taxid['gga'] = 9031
    code2taxid['bta'] = 9913
    code2taxid['cel'] = 6239
    code2taxid['ath'] = 3702
    code2taxid['dre'] = 7955
    code2taxid['dme'] = 7227
    code2taxid['sce'] = 4932

    p = parser.parse_args()

    if not os.path.exists(p.pef) or not os.path.exists(p.info) or p.pn == '':
        parser.print_help()
        sys.exit()

    bin = get_bin()

    # taxid = int(p.taxid)
    # if p.org in code2taxid:
    #     taxid = code2taxid[p.org]
    taxid = code2taxid[p.org]
    createdir(p.pn)
    olog = open(p.pn + '/command.log', 'a')

    if p.reanal == 'select':
        print("Run DEP Selection")
        L = get_comparison(p.pn + '/comparison.txt')
        resultdir = p.pn + '/' + p.pn + '_Result'
        if os.path.exists(resultdir):
            shutil.rmtree(resultdir)
        file_selection(p.pf, p.pef, L, p.pn, olog)
        os.system(
            "{} {}/proteinreport.py -i {} -ia {} -t {}".format(get_python(), bin, p.pn + '/' + p.pn + '_Result', p.pn,
                                                               p.reporttype))
        sys.exit()

    if p.reanal == 'report':
        resultdir = p.pn + '/' + p.pn + '_Result'
        if os.path.exists(resultdir+'/html'):
            shutil.rmtree(resultdir+'/html')
        os.system(
            "{} {}/proteinreport.py -i {} -ia {} -t {}".format(get_python(), bin, p.pn + '/' + p.pn + '_Result', p.pn,
                                                               p.reporttype))
        sys.exit()

    # if p.reanal=='ppi':
    #     # 相互作用网络
    #     print("Run Network Analysis")
    #     L = get_comparison(p.pn + '/comparison.txt')
    #     ppi_search2(L, bin, p.pn, taxid, olog)
    #
    #     print("Run file selection")
    #     resultdir = p.pn + '/' + p.pn + '_Result'
    #     if os.path.exists(resultdir):
    #         shutil.rmtree(resultdir)
    #     file_selection(p.pf, p.pef, L, p.pn, olog)

        ###生成报告

        # os.system(
        #     "{} {}/proteinreport.py -i {} -ia {} -t {}".format(get_python(), bin, p.pn + '/' + p.pn + '_Result', p.pn,
        #                                                        p.reporttype))
        # exit()

    paralog = open(p.pn+'/parameter.txt', 'w')
    paralog.write("{}\n".format(' '.join(sys.argv)))
    paralog.write('FC\t{}\n'.format(p.fc))
    paralog.write('pvalue\t{}\n'.format(p.pvalue))
    paralog.write('org\t{}\n'.format(p.org))
    paralog.write('celllocation\t{}\n'.format(p.celllocation))
    paralog.write('taxid\t{}\n'.format(taxid))
    paralog.close()



    print("\n### Missing Value Analysis ###\n")
    miss_value_handling(p.pef, p.info, p.pn,p.isft)

    print("### Get Probability Files ###")
    olog.write("Rscript {}/getProbability_v2.R {} {} {} {} {}".format(bin, p.info, p.pef, p.pn + "/1.Info/missValue_imputation.xlsx", p.pn, p.isft))
    os.system("Rscript {}/getProbability_v2.R {} {} {} {} {}".format(bin, p.info, p.pef, p.pn + "/1.Info/missValue_imputation.xlsx", p.pn, p.isft))

    print("\n### Get All Motif Sequence ###\n")
    get_all_motif_seq(p.pn + '/1.Info/ProbabilitySite.xlsx', p.pn, p.isft)

    #######运行参数解析
    olog.write("Rscript {0}/readParameterPTM.R {1} {2} {3}\n".format(bin, p.info, p.pn+"/1.Info/matrix.txt", p.pn))
    os.system("Rscript {0}/readParameterPTM.R {1} {2} {3}".format(bin, p.info, p.pn+"/1.Info/matrix.txt", p.pn))

    #######提取输入文件
    olog.write("Rscript {0}/extractExpression2.R {1} {2} {3} {4}\n".format(bin, p.pn + "/1.Info/ModifiedProtein.xlsx", p.info, p.pn, p.org))
    os.system("Rscript {0}/extractExpression2.R {1} {2} {3} {4}".format(bin, p.pn + "/1.Info/ModifiedProtein.xlsx", p.info, p.pn, p.org))

    #######运行整体绘图
    print("\n### Run Sample Analysis ###\n")
    del_min = ''
    if p.maskminvalue is False:
        del_min = 'FALSE'
    else:
        del_min = 'TRUE'
    run_sample_analysis(p.pn+'/groupFile.txt', p.scaleMethod, p.col, p.pn, bin, olog, del_min)

    #######判断是不是多组数据，如果是多组数据需要运行ANOVA检验和绘图, ANOVA检验后对结果进行K-means聚类绘图
    #读取anova行
    groupNumb = get_group_numb(p.pn + '/groupFile.txt')
    ipara = open(p.pn+'/comparison.txt')
    for line in ipara.readlines()[1:]:
        row = line.rstrip().split('\t')
        row[0] = row[0].strip()
        if row[0] == 'anova':
            grouporder = row[1].replace(';',',')
            row[1] = row[1].replace(';', '_')
            anova_input_file = p.pn+'/Anova/'+row[1]+'/input.txt'
            anova_group_file = p.pn+'/Anova/'+row[1]+'/group.txt'
            # ANOVA检验
            sig_number = run_anova_analysis(bin, p.pn+'/Anova/' + row[1], p.org, p.celllocation, anova_input_file, anova_group_file, grouporder)
            ianova = open(p.pn+'/Anova/anova_sig_number.txt', 'a')
            ianova.write("{}\t{}\n".format(row[1], sig_number))
            ianova.close()

    #######全部蛋白分析#######
    print("\n### Run All Proteins Analysis ###\n")
    if p.org in supp and os.path.exists(p.pn+'/input.txt'):
        run_all_protein_analysis(bin, p.pn, p.org, p.celllocation)
    else:
        print("有可能是物种参数设置错误，无法分析")
        sys.exit()

    #######全部Motif分析(ProbALL)#######
    print("\n### Run All Probability MotifAnalysis ###\n")
    of = os.path.join(os.getcwd(), p.pn + '/MotifAnalysis')
    print(of)
    if not os.path.exists(of):
        os.mkdir(of)
    if os.path.exists(of+'/All_MotifInput.txt'):
        codeContent='momo motifx -oc /meme --verbosity 1 --width 13 --eliminate-repeats 13 --min-occurrences 5'
        subprocess.run("echo \"welovebp1188\" |sudo -S docker run -v {workpath}:/meme --user `id -u`:`id -g` {docker_name} {codeC}  /meme/All_MotifInput.txt"
                       .format(workpath=of, docker_name='memesuite/memesuite:5.3.3', codeC=codeContent), shell=True, check=True)
    else:
        print("未检测到Motif分析输入文件！")
        sys.exit()

    #######每一组的差异筛选#######
    print("\n### Run DEP Selection ###\n")
    L = get_comparison(p.pn + '/comparison.txt')
    print(L)

    # 各组差异统计分析 #
    run_dep_selection(L, bin, p.pn, olog, p.fc, p.pvalue)

    #######差异Motif分析#######
    print("\n### Run DEP Motif Analysis ###\n")
    get_DEP_motif(L, p.pn)


    olog.write("Rscript {0}/getAllDEP_Result.R {1} {2} {3}\n".format(bin,p.pf, p.isft, p.pn))
    os.system("Rscript {0}/getAllDEP_Result.R {1} {2} {3}".format(bin, p.pf, p.isft, p.pn))

    # 差异蛋白功能分析 #
    print("\n### Run DEP Functional Analysis ###\n")
    get_DEP_protein_id(L, p.pn, p.type)
    run_dep_function_analysis(L, bin, p.pn, p.org, p.celllocation, supp, olog)
    # 差异柱状图合并 #
    merge_updown_bar(L, p.pn, bin, olog)
    #venn分析
    venn_anal(L, p.pn, bin, olog)

    #######挑选分析结果#######
    print("\n### Run file selection ###\n")
    file_selection(p.pef, L, p.pn, olog)
    #
    # ###生成报告
    # os.system("{} {}/proteinreport.py -i {} -ia {} -t {}".format(get_python(), bin, p.pn + '/' + p.pn + '_Result', p.pn, p.reporttype))

    olog.close()

