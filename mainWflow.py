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
# import check
# import wfmo


### MissValue Analysis
def miss_value_handling(input1,para,pn):
    of = pn + "/1.Info"
    if not os.path.exists(of):
        os.mkdir(of)
    if os.path.exists(input1):
        df = pd.read_excel(input1) # sites.xlsx
        gdf = pd.read_excel(para, sheet_name='Sheet1') # parameter.xlsx

        samples = ["PTM.CollapseKey"] + gdf.loc[:, 'Sample'].tolist()
        print(samples)
        dfsub = df.loc[:, samples]
        dfsub.set_index('PTM.CollapseKey')

        # 将样本分组文件转换为dict
        d_group = {}
        for i in range(len(gdf)):
            samp = gdf.iloc[i, 0]
            group = gdf.iloc[i, 1]
            if group not in d_group:
                d_group[group] = [samp]
            else:
                d_group[group].append(samp)
        print(d_group)
        df_data_copy = dfsub.copy(deep=True)
        df_data_copy = df_data_copy.replace('Filtered', np.NaN)

        def Confidence(df=None, d_group=None, confidence=0.5, how="outer", ):
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
            df_filter.to_csv(of + "/filter.csv", index=0)
            return df_

        df2 = Confidence(df=df_data_copy, d_group=d_group, confidence=0.5, how='outer')
        df2.to_excel(of + "/missValue_clean.xlsx", index=None)

        # ???? 参数控制：KNN填充，最小值填充，均值填充等等
        ## KNN method imputer
        imputer = KNNImputer(n_neighbors=10)
        df2.iloc[:, 1:] = imputer.fit_transform(df2.iloc[:, 1:])

        print(df2.head(5).to_string())
        df2.to_excel(of + "/missValue_imputation.xlsx", index=None)


# Sample Analysis
def run_sample_analysis(groupFile, scaleMethod, col, pn, bin, olog, delmin, exp='DIA'):
    of = pn + "/2.SampleAnalysis"
    if not os.path.exists(of):
        os.mkdir(of)
    if os.path.exists(groupFile):
        # # 抽提数据
        # olog.write("Rscript {0}/extractTable.R {1} {2}\n".format(bin, input, of + "/matrix.txt"))
        # os.system("Rscript {0}/extractTable.R {1} {2}".format(bin, input, of + "/matrix.txt"))
        # 探索性绘图
        print("### Run matrix drawing")
        olog.write(
            "Rscript {0}/matrixDraw.R {1} {2} {3} {4} Intensity log2 {5}\n".format(
                bin, pn + "/1.Info/matrix.txt", groupFile, of, col, delmin))
        os.system(
            "Rscript {0}/matrixDraw.R {1} {2} {3} {4} Intensity log2 {5}".format(
                bin, pn + "/1.Info/matrix.txt", groupFile, of, col, delmin))

        olog.write("Rscript {}/svdiPCA2.R {} {} {} uv\n".format(bin, pn + "/1.Info/matrix.txt", groupFile, of))
        os.system("Rscript {}/svdiPCA2.R {} {} {} uv".format(bin, pn + "/1.Info/matrix.txt", groupFile, of))


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



"""全部蛋白的分析"""
def run_all_protein_analysis(bin, pn, org, celllocation):
    print("### Run All Proteins Analysis ###")
    olog.write("Rscript {}/FuncAnal4ID.R -i {} -o {} -s {} -p {} -c 4\n".format(bin, pn + '/input.txt', pn + '/3.AllProtein', org, bin))
    os.system("Rscript {}/FuncAnal4ID.R -i {} -o {} -s {} -p {} -c 4".format(bin, pn + '/input.txt', pn + '/3.AllProtein', org, bin))

    ########细胞定位的分析的绘图
    if os.path.exists(p.pn + '/3.AllProtein/GOEnrich_Species/CC.txt'):
        olog.write("Rscript {}/dCellLoc.R {} {} {}\n".format(bin, pn + '/3.AllProtein/GOEnrich_Species/CC.txt',
                                                             pn + '/3.AllProtein/GO_subcellular_localization', celllocation))
        os.system("Rscript {}/dCellLoc.R {} {} {}".format(bin, pn + '/3.AllProtein/GOEnrich_Species/CC.txt',
                                                          pn + '/3.AllProtein/GO_subcellular_localization', celllocation))

    #####增加KEGG分类信息
    olog.write(
        "{}/addURLClassKEGG {}/3.AllProtein/KEGG/KEGG.enriched.txt {}/3.AllProtein/KEGG/KEGG.enriched.xls\n".format(bin,pn, pn))
    os.system(
        "{}/addURLClassKEGG {}/3.AllProtein/KEGG/KEGG.enriched.txt {}/3.AllProtein/KEGG/KEGG.enriched.xls".format(bin, pn,pn))

    olog.write(
        "Rscript {}/txt2xlsx.R {}/3.AllProtein/KEGG/KEGG.enriched.xls KEGG {}/3.AllProtein/KEGG/KEGG.enriched.xlsx\n".format(
            bin, pn, pn))
    os.system(
        "Rscript {}/txt2xlsx.R {}/3.AllProtein/KEGG/KEGG.enriched.xls KEGG {}/3.AllProtein/KEGG/KEGG.enriched.xlsx".format(
            bin, pn, pn))

    #####所有蛋白结果的绘图
    olog.write("Rscript {}/PAAdraw.R {}\n".format(bin, pn + '/3.AllProtein'))
    os.system("Rscript {}/PAAdraw.R {}".format(bin, pn + '/3.AllProtein'))

    ############挑选出GO的对应表
    olog.write("Rscript {}/selectGO.R {} 4 {} {}\n".format(bin, pn + '/input.txt', org, pn + "/3.AllProtein/GO"))
    os.system("Rscript {}/selectGO.R {} 4 {} {}".format(bin, pn + '/input.txt', org, pn + "/3.AllProtein/GO"))

    ###########生成Summary文件
    olog.write("{}/makeSummary -i {}\n".format(bin, pn))
    os.system("{}/makeSummary -i {}".format(bin, pn))

    olog.write("Rscript {}/tsv2xlsx.R {} {} {}\n".format(bin, pn + '/protein_summary.txt', 'summary', pn + '/summary.xlsx'))
    os.system("Rscript {}/tsv2xlsx.R {} {} {}".format(bin, pn + '/protein_summary.txt', 'summary', pn + '/summary.xlsx'))


## 未仔细修改
"""多组之间比较的方差分析"""
def run_anova_analysis(bin, pn, org, celllocation, input_file, input_group_file, grouporder):
    olog.write("Rscript {}/extractTable.R {} {}\n".format(bin, input_file, pn + '/matrix.txt'))
    os.system("Rscript {}/extractTable.R {} {}".format(bin, input_file, pn + '/matrix.txt'))

    olog.write("Rscript {}/kmeans.R {} {} {} {}\n".format(bin, input_group_file, pn + '/matrix.txt',
                                                           pn + '/kmeans', grouporder))
    os.system("Rscript {}/kmeans.R {} {} {} {}".format(bin, input_group_file, pn + '/matrix.txt',
                                                        pn + '/kmeans', grouporder))

    # 挑选做比较绘图的文件
    olog.write("Rscript {}/selectFileFromTotal.R {} {} {}\n".format(bin, pn + '/kmeans/ANOVA_result.xlsx',
                                                                    pn + '/input.txt', pn + '/ANOVA.sig.txt'))
    os.system("Rscript {}/selectFileFromTotal.R {} {} {}".format(bin, pn + '/kmeans/ANOVA_result.xlsx',
                                                                 pn + '/input.txt', pn + '/ANOVA.sig.txt'))

    # 比较绘图
    olog.write("Rscript {}/dDEP4multi.R {} {} {}\n".format(bin, pn + '/ANOVA.sig.txt', input_group_file,
                                                           pn + '/Multi_group_expression'))
    os.system("Rscript {}/dDEP4multi.R {} {} {}".format(bin, pn + '/ANOVA.sig.txt', input_group_file,
                                                        pn + '/Multi_group_expression'))

    #热图绘制
    olog.write("Rscript {}/drawHeatmap.R {} {} {} multi\n".format(bin, pn + '/ANOVA.sig.txt',input_group_file, pn))
    os.system("Rscript {}/drawHeatmap.R {} {} {} multi".format(bin, pn + '/ANOVA.sig.txt',input_group_file, pn))

    # 显著的蛋白的功能分析和绘图
    olog.write(
        "Rscript {}/FuncAnal4ID.R -i {} -o {} -s {} -p {} -b {} -d 4 -c 4\n".format(bin, pn + '/ANOVA.sig.txt', pn,
                                                                         org, bin, input_file))
    os.system(
        "Rscript {}/FuncAnal4ID.R -i {} -o {} -s {} -p {} -b {} -d 4 -c 4".format(bin, pn + '/ANOVA.sig.txt', pn,
                                                                       org, bin, input_file))

    ########细胞定位的分析的绘图
    if os.path.exists(p.pn + '/AllProtein/GOEnrich_Species/CC.txt'):
        olog.write("Rscript {}/dCellLoc.R {} {} {}\n".format(bin, pn + '/GOEnrich_Species/CC.txt',
                                                             pn + '/GO_subcellular_localization', celllocation))
        os.system("Rscript {}/dCellLoc.R {} {} {}".format(bin, pn + '/GOEnrich_Species/CC.txt',
                                                          pn + '/GO_subcellular_localization', celllocation))

    #####增加KEGG分类信息
    olog.write(
        "{}/addURLClassKEGG {}/KEGG/KEGG.enriched.txt {}/KEGG/KEGG.enriched.xls\n".format(bin, pn, pn))
    os.system(
        "{}/addURLClassKEGG {}/KEGG/KEGG.enriched.txt {}/KEGG/KEGG.enriched.xls".format(bin, pn, pn))

    olog.write(
        "Rscript {}/txt2xlsx.R {}/KEGG/KEGG.enriched.xls KEGG {}/KEGG/KEGG.enriched.xlsx\n".format(
            bin, pn, pn))
    os.system(
        "Rscript {}/txt2xlsx.R {}/KEGG/KEGG.enriched.xls KEGG {}/KEGG/KEGG.enriched.xlsx".format(
            bin, pn, pn))

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

    """
    if platform.system() == 'Windows':
        bin = get_absdir()
    else:
        bin = '/usr/local/bio/PAA/bin'
    """
    taxid = int(p.taxid)
    if p.org in code2taxid:
        taxid = code2taxid[p.org]

    #createdir(p.pn)

    olog = open(p.pn + '/command.log', 'a')

    # if p.reanal == 'select':
    #     print("Run DEP Selection")
    #     L = get_comparison(p.pn + '/comparison.txt')
    #     resultdir = p.pn + '/' + p.pn + '_Result'
    #     if os.path.exists(resultdir):
    #         shutil.rmtree(resultdir)
    #     file_selection(p.pf, p.pef, L, p.pn, olog)
    #     os.system(
    #         "{} {}/proteinreport.py -i {} -ia {} -t {}".format(get_python(), bin, p.pn + '/' + p.pn + '_Result', p.pn,
    #                                                            p.reporttype))
    #     sys.exit()
    #
    # if p.reanal == 'report':
    #     resultdir = p.pn + '/' + p.pn + '_Result'
    #     if os.path.exists(resultdir+'/html'):
    #         shutil.rmtree(resultdir+'/html')
    #     os.system(
    #         "{} {}/proteinreport.py -i {} -ia {} -t {}".format(get_python(), bin, p.pn + '/' + p.pn + '_Result', p.pn,
    #                                                            p.reporttype))
    #     sys.exit()
    #
    # if p.reanal=='ppi':
    #     # 相互作用网络
    #     print("Run Network Analysis")
    #     # ppi_search(L, bin, p.pn, taxid, olog)
    #     L = get_comparison(p.pn + '/comparison.txt')
    #     ppi_search2(L, bin, p.pn, taxid, olog)
    #
    #     print("Run file selection")
    #     resultdir = p.pn + '/' + p.pn + '_Result'
    #     if os.path.exists(resultdir):
    #         shutil.rmtree(resultdir)
    #     file_selection(p.pf, p.pef, L, p.pn, olog)
    #
    #     ###生成报告
    #
    #     os.system(
    #         "{} {}/proteinreport.py -i {} -ia {} -t {}".format(get_python(), bin, p.pn + '/' + p.pn + '_Result', p.pn,
    #                                                            p.reporttype))
    #     exit()



    paralog = open(p.pn+'/parameter.txt', 'w')
    paralog.write("{}\n".format(' '.join(sys.argv)))
    paralog.write('FC\t{}\n'.format(p.fc))
    paralog.write('pvalue\t{}\n'.format(p.pvalue))
    paralog.write('org\t{}\n'.format(p.org))
    paralog.write('celllocation\t{}\n'.format(p.celllocation))
    paralog.write('taxid\t{}\n'.format(taxid))
    paralog.close()

    print("### Miss Value Handling ###")
    miss_value_handling(p.pef, p.info, p.pn)

    print("### Get Probability Information ###")
    olog.write("Rscript {}/getProbability.R {} {} {} {}".format(bin, p.info, p.pef, p.pn + "/1.Info/missValue_imputation.xlsx", p.pn))
    os.system("Rscript {}/getProbability.R {} {} {} {}".format(bin, p.info, p.pef, p.pn + "/1.Info/missValue_imputation.xlsx", p.pn))

    #######运行参数解析
    olog.write("Rscript {0}/readParameterPTM.R {1} {2} {3}\n".format(bin, p.info, p.pn+"/1.Info/matrix.txt", p.pn))
    os.system("Rscript {0}/readParameterPTM.R {1} {2} {3}".format(bin, p.info, p.pn+"/1.Info/matrix.txt", p.pn))

    #######提取输入文件
    olog.write("Rscript {0}/extractExpression2.R {1} {2} {3} {4}\n".format(bin, p.pn + "/1.Info/ModifiedProtein.xlsx", p.info, p.pn, p.org))
    os.system("Rscript {0}/extractExpression2.R {1} {2} {3} {4}".format(bin, p.pn + "/1.Info/ModifiedProtein.xlsx", p.info, p.pn, p.org))

    #######运行整体绘图
    print("### Run Sample Drawing ###")
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

    #######全部蛋白分析
    if p.org in supp and os.path.exists(p.pn+'/input.txt'):
        run_all_protein_analysis(bin, p.pn, p.org, p.celllocation)
    else:
        print("有可能是物种参数设置错误，无法分析")
        sys.exit()

    #######全部Motif分析(ProbALL)
    of = p.pn + '/MotifAnalysis'
    if not os.path.exists(of):
        os.mkdir(of)
    if os.path.exists(of+'/All_MotifInput.txt'):
        codeContent='momo motifx -oc /meme --verbosity 1 --width 13 --eliminate-repeats 13 --min-occurrences 5'
        subprocess.run("sudo docker run -v {workpath}:/meme --user `id -u`:`id -g` {docker_name} {codeC}  /meme/All_MotifInput.txt"
                       .format(workpath=of, docker_name='memesuite/memesuite:5.3.3', codeC=codeContent), shell=True, check=True)

    else:
        print("未检测到Motif分析输入文件！")
        sys.exit()

