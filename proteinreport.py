#!/usr/bin/python3
import os
import sys
from wfkit import cpfile, get_absdir
import platform
import shutil
import argparse
import hylicense.eda

class protreport:
    def __init__(self, pn, output, resultpn):
        self.pn = pn
        self.resultpn = resultpn
        self.encode = self.get_encode()
        if os.path.exists(self.resultpn + '/multiUpdown.stat.txt'):
            self.comp_group, self.updown_stat = self.load_updown_stat(self.resultpn + '/multiUpdown.stat.txt')
        else:
            print("缺少上下调蛋白质的统计文件，不能完成报告")
            sys.exit()
        self.htmlout = self.pn + '/html'
        if not os.path.exists(self.htmlout):
            os.makedirs(self.htmlout)
        self.output = output
        self.ofile = open(self.output, 'w')
        self.ofile.write('''<html>
<frameset cols="400,*">
<frame src="./html/left.html">
<frame src="./html/firstpage.html" name="rightframe">
</frameset>

</html>
''')
        self.ofile.close()
        self.index = 1

        self.lefthandle = open(self.pn + '/html/left.html','w')
        self.print_head_left(self.lefthandle)
        self.para = self.get_parameter()

    def get_encode(self):
        if platform.system() == 'Windows':
            encode = 'GB2312'
        else:
            encode = 'utf-8'
        return encode

    def get_parameter(self):
        para = dict()
        if os.path.exists(self.resultpn + '/parameter.txt'):
            ipara = open(self.resultpn + '/parameter.txt')
            for line in ipara.readlines():
                row = line.rstrip().split('\t')
                if len(row) > 1:
                    para[row[0]] = row[1]

        return para



    """获取总蛋白数"""
    def get_total_protein_numb(self):
        ifile = open(self.resultpn+'/input.txt')
        totalnumb = len(ifile.readlines()[1:])
        return(totalnumb)

    """获取每组的上下调蛋白数"""
    def load_updown_stat(self,updown_stat_file):
        ifile = open(updown_stat_file)
        comp_group=set()
        udstat = dict()
        for line in ifile.readlines()[1:]:
            row = line.rstrip().split('\t')
            comp_group.add(row[2])
            print(row[2] + '\t' + row[0])
            udstat[row[2] + '\t' + row[0]] = int(row[1])
        return(list(comp_group), udstat)

    def write_image_html(self, image_url, image_caption, image_width=800):
        self.ofile.write('''<div class="text-center">
        <figure class="figure">
                <a href="{0}" class="thumbnail" target=_blank>
                    <img src="{0}" class="figure-img img-fluid d-block mx-auto" width={2} alt="无该结果"  />
                </a>
                <figcaption class="figure-caption text-center">{1}</figcaption>
</figure></div>'''.format(image_url, image_caption, image_width))

    def write_h1_html(self, title):
        self.ofile.write("<h1>{} {}</h1>".format(self.index, title))

    def print_head(self, handle):
        handle.write('''<html><head>
         <meta charset="{}">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <link rel="stylesheet" href="https://cdn.staticfile.org/twitter-bootstrap/4.3.1/css/bootstrap.min.css">
  <link rel="stylesheet" href="CSS/style.css">
  <script src="https://cdn.staticfile.org/jquery/3.2.1/jquery.min.js"></script>
  <script src="https://cdn.staticfile.org/popper.js/1.15.0/umd/popper.min.js"></script>
  <script src="https://cdn.staticfile.org/twitter-bootstrap/4.3.1/js/bootstrap.min.js"></script>
          \n<title>上海拜谱生物蛋白质组学分析报告</title>
        \n</head>
    \n<body><div class="container"><div class="text-right" ><img src="img/baipu.gif" height=50 /></div>'''.format(self.encode))

    def print_head_left(self, handle):
        handle.write('''<html><head>
           <meta charset="{}">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <link rel="stylesheet" href="https://cdn.staticfile.org/twitter-bootstrap/4.3.1/css/bootstrap.min.css">
    <link rel="stylesheet" href="CSS/style2.css">
    <script src="https://cdn.staticfile.org/jquery/3.2.1/jquery.min.js"></script>
    <script src="https://cdn.staticfile.org/popper.js/1.15.0/umd/popper.min.js"></script>
    <script src="https://cdn.staticfile.org/twitter-bootstrap/4.3.1/js/bootstrap.min.js"></script>
            \n<title>上海拜谱生物蛋白质组学分析报告</title>
          \n</head>
      \n<body><div class="container">'''.format(self.encode))

    def print_bottom(self, handle):
        handle.write("""<div class="gohome">回到顶部</div>
  <script src="js/back.js"></script><div class="pageend text-dark">上海拜谱生物科技有限公司 | Shanghai Bioprofile Biotechnology Co., Ltd
Add: 上海市闵行区紫月路468号5楼 | Web: www.bioprofile.cn | Tel: 4008208531 | Email: info@bioprofile.cn</div></div></body></html>""")

    def print_image(self, handle, image_url, image_caption, image_width=600):
        handle.write('''<div class="text-center">
                <figure class="figure">
                        <a href="{0}" class="thumbnail" target=_blank>
                            <img src="{0}" class="figure-img img-fluid d-block mx-auto" width={2} alt="无该结果"  />
                        </a>
                        <figcaption class="figure-caption text-center">{1}</figcaption>
        </figure></div>'''.format(image_url, image_caption, image_width))

    def return_image(self, image_url, image_caption, image_width=600):
        s = '''<div class="text-center">
                <figure class="figure">
                        <a href="{0}" class="thumbnail" target=_blank>
                            <img src="{0}" class="figure-img img-fluid d-block mx-auto" width={2} alt="无该结果"  />
                        </a>
                        <figcaption class="figure-caption text-center">{1}</figcaption>
        </figure></div>'''.format(image_url, image_caption, image_width)
        return s

    def print_h1(self, handle, title, url, name):
        handle.write("<h1>{}. {}<a name=\"{}\"></a></h1>".format(self.index, title,name))
        self.lefthandle.write("""<div class="level1"><a href="{}#{}" target="rightframe">{}</a></div>\n""".format(url, name, str(self.index) + '. ' +title))

    def print_h2(self, handle, title, h2id, url, name):
        handle.write("<h2>{}.{} {}<a name=\"{}\"></a></h2>".format(self.index, h2id, title,name))
        self.lefthandle.write("""<div class="level2"><a href="{}#{}" target="rightframe">{}</a></div>\n""".format(url, name, str(self.index)+"."+str(h2id)+' '+title))

    def print_h3(self, handle,title, h2id, h3id, url, name):
        writetitle = str(self.index) + '.' + str(h2id) + '.' + str(h3id)
        handle.write("<h3>{}.{}.{} {}<a name=\"{}\"></a></h3>".format(self.index, h2id, h3id, title, name))
        self.lefthandle.write("""<div class="level3"><a href="{}#{}" target="rightframe">{}</a></div>\n""".format(url, name, writetitle))

    def print_accordion_head(self,handle, id):
        handle.write("""<div class="accordion" id="{}">""".format(id))

    def print_accordion_bottom(self, handle):
        handle.write("</div>")

    def print_accordion_element(self, handle, accord_id, card_header_id, colapse_id, title, content, show):
        handle.write("""<div class="card">
    <div class="card-header" id="{0}">
      <h2 class="mb-0">
        <button class="btn btn-link btn-block text-left" type="button" data-toggle="collapse" data-target="#{1}" aria-expanded="true" aria-controls="{1}">
          {3}
        </button>
      </h2>
    </div>

    <div id="{1}" class="collapse {5}" aria-labelledby="{0}" data-parent="#{2}">
      <div class="card-body">
        {4}
      </div>
    </div>
  </div>
        """.format(card_header_id, colapse_id, accord_id, title, content,show))

    def print_all_protein_function_analysis(self):
        apo = open(self.htmlout+'/allprot.html','w')
        self.print_head(apo)

        if os.path.exists(self.pn+'/4.AllProtein'):
            apo.write("<div id=\"allprotein\">")
            self.print_h1(apo, """全蛋白的功能分析""", "allprot.html","allprot")
            h2id = 1
            figid = 1

            #apo.write("<h2>{} <a name=\"allprot\"></a>全蛋白GO功能注释与富集（模式生物）</h2>".format(str(self.index)+'.'+str(h2id)))
            self.print_h2(apo, "全蛋白GO功能注释与富集（模式生物）", h2id, "allprot.html", "allprot{}".format(h2id))
            h2id = h2id + 1

            apo.write('''<div><p>GO数据库（The Gene Ontology knowledgebase）是全球最大的有关基因功能的信息来源。在生物医学研究中，GO数据库可作为分子生物学和遗传学实验的大数据计算分析的基础。GO数据库将基因和基因产物从3个方面（Term）进行描述：生物过程（Biological Process，BP）、分子功能（Molecular Function，MF）和细胞组分（Cellular Component，CC）。这3个term就是GO注释的第一层级（level 1），每个分支又会分出更具体的term，依次为level 2、level 3等。</p></div>''')

            #golevel2
            allprotgo=self.pn+'/4.AllProtein/4.1GO'
            urlallprotgo = '../4.AllProtein/4.1GO'
            if os.path.exists(allprotgo+'/GOlevel2CountBar.png'):
                apo.write('''<div><p>对该项目鉴定到的所有蛋白进行GO注释，下图展示了level 2的注释结果，横坐标为level 2的term，左侧纵坐标为count（注释到该term的蛋白数量），右侧纵坐标为percentage（注释到该term的蛋白数量 / 有GO注释的总数量），不同的颜色分别表示3大分支：生物过程（Biological Process，BP）、分子功能（Molecular Function，MF）和细胞组分（Cellular Component，CC）。每个柱子上方标有具体的count。</p><p>结果目录：<span class=\"text-primary\">4. AllProtein\\4.1GO\\GOlevel2CountBar</span></p><div>''')
                self.print_image(apo, urlallprotgo + '/GOlevel2CountBar.png','图{}-{}. level 2的GO注释结果'.format(self.index, figid))
                figid = figid + 1


            if os.path.exists(allprotgo+'/pCountPoint.png'):
                apo.write('''<div><p>进一步筛选统计学意义上更为重要的term，依据Fisher精确检验的算法对注释结果进行功能富集分析（全蛋白的富集分析仅模式生物可做）。结果用p-value表示，小于0.05为功能显著富集，值越小表示功能富集越显著。下图展示了BP、MF和CC 3个分支下富集显著性分别top 10的term。横坐标为富集显著性p-value的负对数转化，纵坐标为GO term，每个圆圈表示一个term，圆圈大小表示count，不同的颜色分别表示3大分支。</p><p>结果目录：<span class="text-primary">4. AllProtein\\4.1GO\\pCountPoint</span></p></div>''')
                self.print_image(apo,urlallprotgo + '/pCountPoint.png', '图{}-{}. GO功能富集气泡图（count）'.format(self.index, figid))
                figid = figid + 1


            if os.path.exists(allprotgo+'/pRFPoint.png'):
                apo.write('''<div><p>下图也展示了显著性top 10的三大分支term，只是圆圈大小表示富集因子（Rich factor），Rich factor=（a/b）/（c/d），a表示注释到该term的蛋白数，b表示注释到term的蛋白总数，c表示该term的背景蛋白数，d表示总term的背景蛋白数。</p><p>结果目录：<span class="text-primary">4. AllProtein\\4.1GO\\pRFPoint</span></p></div>''')
                self.print_image(apo,urlallprotgo + '/pRFPoint.png', '图{}-{}. GO功能富集气泡图（rich  factor）'.format(self.index, figid))
                figid =figid + 1

            if os.path.exists(allprotgo+'/BP.EnrichedBar.png') or os.path.exists(allprotgo+'/MF.EnrichedBar.png') or os.path.exists(allprotgo+'/CC.EnrichedBar.png'):
                apo.write('''<div><p>接下来对3大分支的富集结果单独做统计。</p></div>''')
                apo.write('''
                <ul class="nav nav-tabs" role="tablist">
                    <li class="nav-item">
                        <a class="nav-link active" data-toggle="tab" href="#bp">Biological Process</a>
                    </li>
                    <li class="nav-item">
                        <a class="nav-link" data-toggle="tab" href="#cc">Cellular Component</a>
                    </li>
                    <li class="nav-item">
                        <a class="nav-link" data-toggle="tab" href="#mf">Molecular Function</a>
                    </li>
                </ul>
                <div class="tab-content">
                    ''')

                #生物过程(Biological Process)
                apo.write('''<div id="bp" class="container tab-pane active"><br>
                <p>下图展示了富集显著性top 10的生物过程(Biological Process,BP)模体富集柱状图，横坐标为富集显著性p-value的负对数转化，纵坐标为BP分支下的term，柱子越长表示该term的富集显著性越高，每个柱子右侧都有具体的count和p-value。</p><p>结果目录：<span class="text-primary">4. AllProtein\\4.1GO\\BP.EnrichedBar</span></p>
                ''')
                self.print_image(apo,urlallprotgo+'/BP.EnrichedBar.png', "图{}-{}. GO的BP功能富集柱状图".format(self.index, figid))
                figid = figid + 1

                apo.write('''<p>下图展示了富集显著性top 10的生物过程(Biological Process,BP)模体富集气泡图。横坐标表示rich factor，纵坐标表示富集显著性p-value的负对数转化。图中每个圆圈表示一个BP term，用不同的颜色表示。圆圈大小表示count。</p><p>结果目录：<span class="text-primary">4. AllProtein\\4.1GO\\BP.RichFactor</span></p>''')
                self.print_image(apo,urlallprotgo+'/BP.RichFactor.png', "图{}-{}. BP term的显著性top 10富集气泡图".format(self.index, figid))
                figid = figid + 1

                apo.write('''<p>下图展示了显著性top20的BP term，同样用气泡图展示。这里用横坐标表示p-value的负对数转化，圆圈颜色表示rich factor，圆圈大小表示count。</p><p>结果目录：<span class="text-primary">4. AllProtein\\4.1GO\\BP.EnrichedSymbol</span></p>''')
                self.print_image(apo,urlallprotgo +'/BP.EnrichedSymbol.png', '图{}-{}. BP term的显著性top 20富集气泡图'.format(self.index, figid))
                figid = figid + 1

                apo.write('''<p>下图为BP功能富集的有向无循环图，对显著性最高的前10个term作为有向无循环图的主节点，以方框表示。并通过包含关系（is a）将相关联的term一起展示，用圆框表示。包含关系（is a）的具体解释可参见http://geneontology.org/docs/ontology-relations/。方框或圆框中展示了term 的ID、description、p值、count与背景数。颜色越接近红色富集显著性越高。</p><p>结果目录：<span class="text-primary">4. AllProtein\\4.1GO\\BP.DAG</span></p>''')
                self.print_image(apo,urlallprotgo+'/BP.DAG.png', "图{}-{}. GO的BP功能富集有向无循环图".format(self.index, figid))
                figid = figid + 1

                apo.write('</div>')

                # Cellular Component
                apo.write('''<div id="cc" class="container tab-pane fade"><br>
                                <p>下图展示了富集显著性top 10的细胞组份(Cellular Component)模体富集柱状图，横坐标为富集显著性p-value的负对数转化，纵坐标为CC分支下的term，柱子越长表示该term的富集显著性越高，每个柱子右侧都有具体的count和p-value。</p>
                                <p>结果目录：<span class="text-primary">4. AllProtein\\4.1GO\\CC.EnrichedBar</span></p>
                                ''')
                self.print_image(apo, urlallprotgo + '/CC.EnrichedBar.png',
                                      "图{}-{}. GO的CC功能富集柱状图".format(self.index, figid))
                figid = figid + 1

                apo.write('''<p>下图展示了富集显著性top 10的细胞组份(Cellular Component)模体富集气泡图。横坐标表示rich factor，纵坐标表示富集显著性p-value的负对数转化。图中每个圆圈表示一个CC term，用不同的颜色表示。圆圈大小表示count。</p>
                                <p>结果目录：<span class="text-primary">4. AllProtein\\4.1GO\\CC.RichFactor</span></p>''')
                self.print_image(apo, urlallprotgo + '/CC.RichFactor.png',
                                      "图{}-{}. CC term的显著性top 10富集气泡图".format(self.index, figid))
                figid = figid + 1

                apo.write('''<p>下图展示了显著性top20的CC term，同样用气泡图展示。这里用横坐标表示p-value的负对数转化，圆圈颜色表示rich factor，圆圈大小表示count。</p>
                                <p>结果目录：<span class="text-primary">4. AllProtein\\4.1GO\\CC.EnrichedSymbol</span></p>
                                ''')
                self.print_image(apo, urlallprotgo + '/CC.EnrichedSymbol.png',
                                      '图{}-{}. CC term的显著性top 20富集气泡图'.format(self.index, figid))
                figid = figid + 1

                apo.write('''<p>下图为CC功能富集的有向无循环图，对显著性最高的前10个term作为有向无循环图的主节点，以方框表示。并通过包含关系（is a）将相关联的term一起展示，用圆框表示。包含关系（is a）的具体解释可参见http://geneontology.org/docs/ontology-relations/。方框或圆框中展示了term 的ID、description、p值、count与背景数。颜色越接近红色富集显著性越高。</p>
                                <p>结果目录：<span class="text-primary">4. AllProtein\\4.1GO\\CC.DAG</span></p>''')
                self.print_image(apo,urlallprotgo + '/CC.DAG.png', "图{}-{}. GO的CC功能富集有向无循环图".format(self.index, figid))
                figid = figid + 1
                apo.write('</div>')

                # Molecular Function
                apo.write('''<div id="mf" class="container tab-pane fade"><br>
                <p>下图展示了富集显著性top 10的分子功能(Molecular Function, MF)模体富集柱状图，横坐标为富集显著性p-value的负对数转化，纵坐标为MF分支下的term，柱子越长表示该term的富集显著性越高，每个柱子右侧都有具体的count和p-value。</p>
                                                <p>结果目录：<span class="text-primary">4. AllProtein\\4.1GO\\MF.EnrichedBar</span></p>
                                                ''')
                self.print_image(apo,urlallprotgo + '/MF.EnrichedBar.png',
                                      "图{}-{}. GO的MF功能富集柱状图".format(self.index, figid))
                figid = figid + 1

                apo.write('''<p>下图展示了富集显著性top 10的分子功能(Molecular Function, MF)模体富集气泡图。横坐标表示rich factor，纵坐标表示富集显著性p-value的负对数转化。
                                                图中每个圆圈表示一个MF term，用不同的颜色表示。圆圈大小表示count。</p>
                                                <p>结果目录：<span class="text-primary">4. AllProtein\\4.1GO\\MF.RichFactor</span></p>''')
                self.print_image(apo,urlallprotgo + '/MF.RichFactor.png',
                                      "图{}-{}. MF term的显著性top 10富集气泡图".format(self.index, figid))
                figid = figid + 1

                apo.write('''<p>下图展示了显著性top20的MF term，同样用气泡图展示。这里用横坐标表示p-value的负对数转化，圆圈颜色表示rich factor，圆圈大小表示count。</p>
                                                <p>结果目录：<span class="text-primary">4. AllProtein\\4.1GO\\MF.EnrichedSymbol</span></p>
                                                ''')
                self.print_image(apo,urlallprotgo + '/MF.EnrichedSymbol.png',
                                      '图{}-{}. MF term的显著性top 20富集气泡图'.format(self.index, figid))
                figid = figid + 1

                apo.write('''<p>下图为MF功能富集的有向无循环图，对显著性最高的前10个term作为有向无循环图的主节点，以方框表示。
                                                并通过包含关系（is a）将相关联的term一起展示，用圆框表示。包含关系（is a）的具体解
                                                释可参见http://geneontology.org/docs/ontology-relations/。
                                                方框或圆框中展示了term 的ID、description、p值、count与背景数。颜色越接近红色富集显著性越高。</p>
                                                <p>结果目录：<span class="text-primary">4. AllProtein\\4.1GO\\CC.DAG</span></p>''')
                self.print_image(apo,urlallprotgo + '/MF.DAG.png', "图{}-{}. GO的CC功能富集有向无循环图".format(self.index, figid))
                figid = figid + 1
                apo.write('</div>')
                apo.write('</div>')

            #全蛋白KEGG结果
            #apo.write("<h2>{} 全蛋白KEGG pathway注释与富集</h2>".format(str(self.index) + '.' + str(h2id)))
            self.print_h2(apo, "全蛋白KEGG pathway注释与富集", h2id, "allprot.html", "allprot{}".format(h2id))

            h2id = h2id + 1
            apo.write('''<div><p>KEGG数据库（Kyoto Encyclopedia of Genes and Genomes）用于从分子水平的信息中了解细胞、生
            物体和生态系统的高级功能和用途，尤其是从基因组测序和其他高通量实验技术生成的大规模分子数据集中获取信息。在生物体中，蛋
            白质并不独立行使其功能，而是不同蛋白质相互协调完成一系列生化反应以行使其生物学功能。因此，KEGG通路（KEGG pathway）分
            析是更系统、全面地了解细胞的生物学过程、性状或疾病的发生机理、药物作用机制等最直接和必要的途径，将有助于我们去了解蛋白的
            功能和蛋白间的关系。</p></div>''')

            allprotkegg = self.pn + '/4.AllProtein/4.2KEGG'
            urlallprokegg = '4.AllProtein/4.2KEGG'
            if os.path.exists(allprotkegg + '/Top10EnrichedBar.png'):
                apo.write('''<div><p>依据Fisher精确检验的算法对注释结果进行功能富集分析（全蛋白的富集分析仅模式生物可做）。结果用
            p-value表示，小于0.05为功能显著富集，值越小表示功能富集越显著。下图展示了富集显著性top10的pathway，横坐标为p-value的
            负对数转化，纵坐标为pathway，柱子右侧标有具体的count（注释到该pathway的蛋白数量）和p-value数值。</p>
            <p>结果目录：<span class="text-primary">4. AllProtein\\4.2KEGG\\Top10EnrichedBar</span></p></div>''')
                self.print_image(apo,"../"+urlallprokegg + '/Top10EnrichedBar.png', "图{}-{}. 全蛋白的KEGG通路富集柱状图".format(self.index, figid))
                figid = figid + 1

            if os.path.exists(allprotkegg + '/Top10EnrichSymbol.png'):
                apo.write('''<div><p>下图展示了富集显著性top 10的pathway，用气泡图展示。横坐标表示功能富集因子
                （rich factor），纵坐标表示富集显著性p-value的负对数转化，图中每个圆圈表示一个pathway，用不同的颜色表示，
                圆圈大小表示count。Rich factor=（a/b）/（c/d），a表示注释到该pathway的蛋白数，b表示注释到pathway的蛋白
                总数，c表示该pathway的背景蛋白数，d表示总pathway的背景蛋白数。</p>
                <p>结果目录：<span class="text-primary">4. AllProtein\\4.2KEGG\\Top10EnrichSymbol</span></p><div>
                ''')
                self.print_image(apo,"../"+urlallprokegg + '/Top10EnrichSymbol.png', "图{}-{}. 全蛋白的KEGG通路富集气泡图（top10）".format(self.index, figid))
                figid = figid + 1

            if os.path.exists(allprotkegg + '/Top20EnrichedSymbol2.png'):
                apo.write('''<div><p>下图展示了显著性top20的pathway，同样用气泡图展示，可以将p-value、rich factor、
                count同时展示在一张图上。这里用横坐标表示p-value的负对数转化，圆圈颜色表示rich factor，圆圈大小表示count。</p>
                <p>结果目录：<span class="text-primary">4. AllProtein\\4.2KEGG\\Top20EnrichedSymbol2</span></p></div>
                ''')
                self.print_image(apo,"../"+urlallprokegg + '/Top20EnrichedSymbol2.png', "图{}-{}. 全蛋白的KEGG通路富集气泡图（top20）".format(self.index, figid))
                figid = figid + 1

            if os.path.exists(allprotkegg + '/KEGG.Sig.Bar.png'):
                apo.write('''<div><p>下图展示了富集显著性top30的KEGG pathway。横坐标为p-value的负对数转化，纵坐标为
                具体通路。柱子颜色归属level 1分类。柱子右侧有具体的count和p-value。</p>
                <p>结果目录：<span class="text-primary">4.AllProtein\\4.2KEGG\\KEGG.Sig.Bar</span></p></div>
                ''')
                self.print_image(apo,"../"+urlallprokegg + '/KEGG.Sig.Bar.png', '''图{}-{}. 全蛋白的KEGG通路富集柱状图<br>
                （不同颜色归属level 1分类的7大分支：代谢（Metabolism，M），遗传信息处理（Genetic Information Processing，G），
                环境信息处理（Environmental Information Processing，E），细胞过程（Cellular Processes，C），
                生物体系统（Organismal Systems，O），人类疾病（Human Diseases，H），药物开发（Drug Development，D））'''.format(self.index, figid))
                figid = figid + 1

            if os.path.exists(allprotkegg + '/KEGG.Sig.Bar2.png'):
                apo.write('''<div><p>下图也展示了富集显著性top30的KEGG pathway。柱子颜色归属level 2分类。</p>
                <p>结果目录：<span class="text-primary">4. AllProtein\\4.2KEGG\\KEGG.Sig.Bar2</span></p></div>
                ''')
                self.print_image(apo,"../"+urlallprokegg + '/KEGG.Sig.Bar2.png', "图{}-{}. 全蛋白的KEGG通路富集柱状图（不同颜色归属level 2 分类）".format(self.index, figid))

            if os.path.exists(allprotkegg + '/kegg.sig.png'):
                apo.write('''<div><p>下图也展示了富集显著性top 30的通路，用气泡图展示。将pathway归属到level 1分类，每个类
                的pathway从上至下-log10（p-value）依次降低，即p-value依次升高，显著性依次降低。圆圈大小表示count，
                颜色表示rich factor。</p>
                <p>结果目录：<span class="text-primary">4. AllProtein\\4.2KEGG\\kegg.sig</span></p></div>
                ''')
                self.print_image(apo,"../"+urlallprokegg + '/kegg.sig.png', "图{}-{}. 全蛋白的KEGG通路富集气泡图".format(self.index, figid))

            #apo.write("<h2>{} 全蛋白亚细胞定位分析（CC）</h2>".format(str(self.index) + '.' + str(h2id)))
            self.print_h2(apo, "全蛋白亚细胞定位分析（CC）", h2id, "allprot.html", "allprot{}".format(h2id))

            h2id = h2id + 1

            apo.write('''<div><p>细胞可以分为多个细胞器或者细胞区域，如细胞核（nucleolus），细胞质（cytoplasm），
            线粒体（mitochondrion），叶绿体（chloroplast），胞外区（extracellular region），核膜（nuclear membrane），
            过氧化物酶体（peroxisome），内质网（endoplasmic reticulum），高尔基体（Golgi apparatus），溶酶体（lysosome），
            液泡（vacuole）等，这些细胞器被称为亚细胞。不同的基因产物只有在特定的亚细胞中才能正常发挥功能，因此预测蛋白亚细胞定位
            对了解蛋白结构和特性以及蛋白间的相互作用具有重要意义。</p></div>''')

            if os.path.exists(self.pn + '/4.AllProtein/4.3subcellular_localization/CellLocationMerge.png'):
                apo.write('''<div><p>这里对GO数据库中的细胞组分（Cellular Component，CC，亚细胞定位）进行注释和统计，下图展示了亚细胞定位的蛋白数量统
                计结果，依次为柱形图、气泡图、饼图和环图。</p>
                <p>结果目录：<span class="text-primary">4. AllProtein\\4.4subcellular localization</span></p></div>''')
                self.print_image(apo,"../"+'4.AllProtein/4.3subcellular_localization/CellLocationMerge.png', "图{}-{}. GO的CC亚细胞定位结果".format(self.index, figid))

            self.index = self.index + 1
            self.print_bottom(apo)
            apo.close()

    def print_dep_protein_function_analysis(self):
        dpo = open(self.htmlout+'/dep.html','w')
        depstat = self.pn+'/5.DEPStatistics'
        depstaturl='../5.DEPStatistics'
        self.print_head(dpo)
        self.print_h1(dpo, "差异蛋白的表达与功能分析", "dep.html", "dep")
        h2id = 1
        figid = 1
        tableid = 1
        #dpo.write("<h2>{}.{} 差异蛋白筛选</h2>".format(self.index, h2id))
        self.print_h2(dpo, "差异蛋白筛选", h2id, "dep.html",
                      "dep{}".format(h2id))

        h2id = h2id + 1
        dpo.write("""<div><p>在蛋白质的显著性差异分析中，采用T检验（Student t test）
        结合变化倍数（Fold change，FC，两组间表达量平均值的比值）的方法，筛选出显著差异蛋白（通常满足p-value<{2}，同时FC>{3}或<{4}）。</p>
        <p>结果目录：<span class="text-primary">5.DEPStatistics\\DEP</span>，其中红色表示显著上调差异蛋白，蓝色表示显著下调差异蛋白。
        统计结果请见表{0}-{1}。</div></p>
        """.format(self.index, tableid, self.para["pvalue"], self.para["FC"], str(1)+'/'+str(self.para["FC"])))

        ''' 上下调统计表'''
        dpo.write("<div>表{}-{}. 差异蛋白质结果统计</div>".format(self.index, tableid))
        tableid = tableid + 1
        totalnumb = self.get_total_protein_numb()
        dpo.write("""<table class="table"><thead><tr>
      <th scope="col">差异比较组</th>
      <th scope="col">上调<br />up-regulation</th>
      <th scope="col">下调down-regulation</th>
      <th scope="col">显著差异蛋白总数</th></tr>
        </thead><tbody>
                    """)
        for comp in sorted(self.comp_group):
            print(comp)
            if comp+'\tUp' in self.updown_stat:
                up_numb = self.updown_stat[comp+'\tUp']
            else:
                up_numb = 0

            if comp+'\tDown' in self.updown_stat:
                down_numb = self.updown_stat[comp+'\tDown']
            else:
                down_numb = 0
            dep_numb = up_numb + down_numb
            dpo.write("""<tr><td>{}</td>
            <td>{}</td>
            <td>{}</td>
            <td>{}</td>
            </tr>
            """.format(comp, up_numb,down_numb,dep_numb))

        dpo.write("</tbody></table>")

        if os.path.exists(depstat+'/multiDEPStatbar.png'):
            self.print_image(dpo, depstaturl+"/multiDEPStatbar.png", "图{}-{} 差异蛋白质数统计的柱状图".format(self.index, figid))
            figid = figid + 1

        #如果存在多组比较，呈现Venn分析
        if os.path.exists(depstat+'/Venn/all.eulerEllipse.png'):
            #dpo.write("<h2>{}.{} Venn分析</h2>".format(self.index, h2id))
            self.print_h2(dpo, "Venn分析", h2id, "dep.html",
                          "dep{}".format(h2id))
            h2id = h2id + 1
            dpo.write("""<div><p>采用Venn分析展示不同比较组差异筛选得到的蛋白之间的交叠情况。</p>
            <p>结果目录：<span class="text-primary">5.DEPStatistics\\Venn<span></p></div>""")
            self.print_image(dpo, depstaturl+'/Venn/all.eulerEllipse.png', "图{}-{} 比较组筛选到差异蛋白的维恩图".format(self.index, figid))
            figid = figid + 1

        #火山图
        """组别数"""
        comp_numb = len(self.comp_group)
        #用折叠显示火山图
        #dpo.write("<h2>{}.{} 火山图</h2>".format(self.index, h2id))
        self.print_h2(dpo, "火山图", h2id, "dep.html",
                      "dep{}".format(h2id))
        h2id = h2id + 1
        dpo.write("""\n<div><p>横坐标为FC的对数转化，纵坐标为p-value的负对数转化，每个点表示一个蛋白，其中红色点为显著高表达蛋白质，
                  颜色越深上调倍数越高；蓝色点为显著低表达蛋白质，颜色越深下调倍数越高；灰色点为非差异蛋白质。</p>
                  <p>结果目录：<span class="text-primary">5. DEPStatistics\\*.vs.*.Volcano</span></p>
    </div>""")
        self.print_accordion_head(dpo, "volcano")
        c=0
        for comp in self.comp_group:
            if c == 0:
                show = 'show'
            else:
                show = ''
            self.print_accordion_element(dpo, "volcano", "volcano" + str(c),
                                         "vol"+str(c), comp, self.return_image(depstaturl+'/'+comp+'/'+comp+'.volcano.png', "图{}-{} {}组的火山图".format(self.index, figid, comp)),show
                                         )
            figid = figid + 1
            c = c+1

        self.print_accordion_bottom(dpo)

        #dpo.write("<h2>{}.{} 差异蛋白聚类分析</h2>".format(self.index, h2id))
        self.print_h2(dpo, "差异蛋白聚类分析", h2id, "dep.html",
                      "dep{}".format(h2id))
        h2id = h2id + 1
        dpo.write("""<div><p>聚类分析是一种常用的探索性数据分析方法，其目的是在相似性的基础上对数据进行分组、归类。聚类树枝越近表示表达
        模式相似性较高，聚类树枝越远表示表达模式相似性较低。聚类算法会对样本（Sample）和变量（Variable，此指蛋白质的定量信息）两个维度进
        行分类。对样本的聚类结果可以检验所筛选的目标蛋白质在组内样本中的表达稳定性；对目标蛋白质的聚类，相近表达模式的蛋白质可能具有相似的
        功能或者参与相同的生物学途径，或者在通路中处于临近的调控位置。</p>
        <p>下图展示了差异蛋白在两组样本中的聚类热图，颜色越红表示相对表达量越高，越蓝表示相对表达量越低。</p>
        <p>基因名称的热图结果目录：<span class="text-primary">5. DEPStatistics\\*.vs.*.genesymbol.pheatmap.png</span></p>
        <p>UniProt蛋白号的热图结果目录：<span class="text-primary">5. DEPStatistics\\*.vs.*.accession.pheatmap.png</span></p>
        </div>
        """)
        self.print_accordion_head(dpo, "heatmap")
        c = 0
        for comp in self.comp_group:
            if c == 0:
                show = 'show'
            else:
                show = ''

            if os.path.exists(depstat + '/' + comp + '/' + comp + '.genesymbol.pheatmap.png'):
                content = self.return_image(depstaturl + '/' + comp + '/' + comp + '.genesymbol.pheatmap.png',
                                            "图{}-{} {}组的热图（基因名称）".format(self.index, figid, comp))
                figid = figid + 1
            else:
                content = ''
            content = content + self.return_image(depstaturl + '/' + comp + '/' + comp + '.accession.pheatmap.png',
                                                           "图{}-{} {}组的热图（UniProt号）".format(self.index, figid, comp))

            self.print_accordion_element(dpo, "heatmap", "heatmap" + str(c),
                                         "heat" + str(c), comp, content, show
                                         )
            figid = figid + 1
            c = c + 1

        self.print_accordion_bottom(dpo)

        #差异蛋白丰度分析
        #dpo.write("<h2>{}.{} 差异蛋白丰度分析</h2>".format(self.index, h2id))
        self.print_h2(dpo, "差异蛋白丰度分析", h2id, "dep.html",
                      "dep{}".format(h2id))
        h2id = h2id + 1
        dpo.write("""<div><p>下图展示了某个差异蛋白在两组中的平均表达情况，横坐标为两个组别，纵坐标为相对表达量。结果提供了柱状图、箱线图
        、抖图和散点图。其中BarComp图中显示的p-value是采用Welch T检验，图中标有具体数值。同时用*表示，*表示p&lt;0.05，**表示p&lt;0.01，
        ***表示p&lt;0.001，ns表示p&gt;>0.05（即没有显著性）。</p>
        <p>结果目录：<span class="text-primary">5. DEPStatistics\\*.vs.*\\*.vs.*.bar、box、jitter、scatter
        </span></p></div>""")
        self.print_image(dpo, depstaturl + '/' + self.comp_group[0] + '/'+ self.comp_group[0] +'.bar/1.png',
                         "图{}-{} 差异蛋白两组表达柱状图".format(self.index, figid),400)
        figid = figid + 1

        self.print_image(dpo, depstaturl + '/' + self.comp_group[0] + '/' + self.comp_group[0] + '.box/1.png',
                         "图{}-{} 差异蛋白两组表达箱线图".format(self.index, figid), 400)
        figid = figid + 1

        self.print_image(dpo, depstaturl + '/' + self.comp_group[0] + '/' + self.comp_group[0] + '.jitter/1.png',
                         "图{}-{} 差异蛋白两组表达抖图".format(self.index, figid), 400)
        figid = figid + 1

        self.print_image(dpo, depstaturl + '/' + self.comp_group[0] + '/' + self.comp_group[0] + '.scatter/1.png',
                         "图{}-{} 差异蛋白两组表达散点图".format(self.index, figid), 400)
        figid = figid + 1

        depfuncurl = '../6.DEPFunction'
        depfun     = self.pn+'/6.DEPFunction'
        #dpo.write("<h2>{}.{} 差异蛋白的GO功能注释与富集</h2>".format(self.index, h2id))
        self.print_h2(dpo, "差异蛋白的GO功能注释与富集", h2id, "dep.html",
                      "dep{}".format(h2id))
        h2id = h2id + 1
        dpo.write("""<div><p>GO数据库（The Gene Ontology knowledgebase）是全球最大的有关基因功能的信息来源。在生物医学研究中，
        GO数据库可作为分子生物学和遗传学实验的大数据计算分析的基础。GO数据库将基因和基因产物从3个方面（Term）进行描述：生物过程
        （Biological Process，BP）、分子功能（Molecular Function，MF）和细胞组分（Cellular Component，CC）。这3个term就是GO
        注释的第一层级（level 1），每个分支又会分出更具体的term，依次为level 2、level 3等。</p>
        <p>对比较组差异蛋白（all diff/ up / down）分别进行GO注释。由于篇幅限制，这里仅展示all diff的统计结果。下图展示了level 2的
        注释结果，横坐标为level 2的term，左侧纵坐标为count（注释到该term的差异蛋白数量），右侧纵坐标为percentage（注释到该term的差
        异蛋白数量 / 有GO注释的差异蛋白总数量），3大分支：生物过程（Biological Process，BP）、分子功能（Molecular
         Function，MF）和细胞组分（Cellular Component，CC）分别用不同的颜色表示。每个柱子上方标有具体的count。</p>
         <p>结果目录：<span class="text-primary">6. DEPFunction\*.vs.*.all_GO\GOlevel2CountBar</span></p></div>
        """)
        self.print_accordion_head(dpo, "golevel")
        c = 0
        for comp in self.comp_group:
            if c == 0:
                show = 'show'
            else:
                show = ''
            self.print_accordion_element(dpo, "golevel", "golevel" + str(c),
                                         "gol" + str(c), comp,
                                         self.return_image(depfuncurl + '/' + comp + '/' + comp + '.all_GO_Level2/GOlevel2CountBar.png',
                                                           "图{}-{} {}组的差异蛋白的GO term注释".format(self.index, figid, comp)), show
                                         )
            figid = figid + 1
            c = c + 1

        self.print_accordion_bottom(dpo)

        #富集分析
        #pCountPoint
        dpo.write("""<div><br /><br /><p>进一步筛选统计学意义上更为重要的term。富集方法为fisher精确检验。结果用p-value表示，小于0.05为功
        能显著富集，值越小表示功能富集越显著。下图展示了BP、MF和CC 3个分支下差异蛋白富集显著性分别top 10的term。横坐标为富集显著性
        p-value的负对数转化，纵坐标为GO term，每个圆圈表示一个term，圆圈大小表示差异蛋白count，3大分支分别用不同的颜色分别表示。我们分别用
        全部鉴定的蛋白质和该物种全部基因为背景参照作了GO功能的富集分析,以下结果均以全部鉴定的蛋白质为背景参照的结果为例</p>
        <p>以全部鉴定的蛋白质为背景参照结果目录：<span class="text-primary">6. DEPFunction\*.vs.*.all_GO_Ident\pCountPoint</span></p>
        <p>以该物种全部基因为背景参照结果目录:<span class="text-primary">6. DEPFunction\*.vs.*.all_GO_Species\pCountPoint</span></p>
        </div>
        """)
        self.print_accordion_head(dpo, "goenrich")
        c = 0
        for comp in self.comp_group:
            if c == 0:
                show = 'show'
            else:
                show = ''
            self.print_accordion_element(dpo, "goenrich", "goenrich" + str(c),
                                         "goe" + str(c), comp,
                                         self.return_image(
                                             depfuncurl + '/' + comp + '/' + comp + '.all_GO_Ident/pCountPoint.png',
                                             "图{}-{} {}组的差异蛋白富集气泡图".format(self.index, figid, comp)), show
                                         )
            figid = figid + 1
            c = c + 1

        self.print_accordion_bottom(dpo)

        #pRFCount
        dpo.write("""<div><br /><br /><p>下图也展示了显著性top 10的三大分支term，只是圆圈大小表示富集因子（Rich factor），Rich factor=（a/b）/（c/d），
        a表示注释到该term的差异蛋白数，b表示注释到term的差异蛋白总数，c表示该term的背景蛋白数，d表示总term的背景蛋白数（模式生物的背景
        蛋白数可以是本项目鉴定到的蛋白，也可以是该物种的蛋白，两个结果都会提供。非模式生物的背景蛋白数只能是本项目鉴定到的蛋白）。</p>
        <p>以全部鉴定的蛋白质为背景参照结果目录：<span class="text-primary">6. DEPFunction\*.vs.*.all_GO_Ident\pRFPoint</span></p>
        <p>以该物种全部基因为背景参照结果目录:<span class="text-primary">6. DEPFunction\*.vs.*.all_GO_Species\pRFPoint</span></p>
        </div>
        """)
        self.print_accordion_head(dpo, "goenrichrf")
        c = 0
        for comp in self.comp_group:
            if c == 0:
                show = 'show'
            else:
                show = ''
            self.print_accordion_element(dpo, "goenrichrf", "goenrichrf" + str(c),
                                         "goerf" + str(c), comp,
                                         self.return_image(
                                             depfuncurl + '/' + comp + '/' + comp + '.all_GO_Ident/pRFPoint.png',
                                             "图{}-{} {}组的差异蛋白富集气泡图（BP，CC and MF）".format(self.index, figid, comp)), show
                                         )
            figid = figid + 1
            c = c + 1

        self.print_accordion_bottom(dpo)

        #updown pvalueBar
        dpo.write("""<div><br /><br /><p>三大分支的富集结果单独进行统计。下图展示了差异蛋白富集显著性top 20的GO term。横坐标为富集显著性p-value的
        负lg转化，纵坐标为term。图中还展示了同一term中的上调&下调蛋白数量比例，用不同颜色表示，每个柱子右侧有标注具体的p值和count。
        <p>以全部鉴定的蛋白质为背景参照结果目录：<span class="text-primary">6. DEPFunction\*.vs.*.all_GO_Ident\[BP|CC|MF].updown.pvalueBar</span></p>
        <p>以该物种全部基因为背景参照结果目录:<span class="text-primary">6. DEPFunction\*.vs.*.all_GO_Species\[BP|CC|MF].updown.pvalueBar</span></p>
        </div>""")
        self.print_accordion_head(dpo, "goupdownbar")
        c = 0
        for comp in self.comp_group:
            if c == 0:
                show = 'show'
            else:
                show = ''
            content=self.return_image(
                                             depfuncurl + '/' + comp + '/' + comp + '.all_GO_Ident/BP.updown.pvalueBar.png',
                                             "图{}-{} {}组的差异蛋白的BP term功能富集bar图（p-value）".format(self.index, figid, comp))
            figid = figid + 1

            content = content + self.return_image(
                                             depfuncurl + '/' + comp + '/' + comp + '.all_GO_Ident/CC.updown.pvalueBar.png',
                                             "图{}-{} {}组的差异蛋白的CC term功能富集bar图（p-value）".format(self.index, figid, comp))
            figid = figid + 1

            content = content + self.return_image(
                                             depfuncurl + '/' + comp + '/' + comp + '.all_GO_Ident/MF.updown.pvalueBar.png',
                                             "图{}-{} {}组的差异蛋白的MF term功能富集bar图（p-value）".format(self.index, figid, comp))
            figid = figid + 1

            self.print_accordion_element(dpo, "goupdownbar", "goupdownbar" + str(c),
                                         "goupd" + str(c), comp,
                                         content, show
                                         )
            figid = figid + 1
            c = c + 1

        self.print_accordion_bottom(dpo)

        #rf
        dpo.write("""<div><br /><br /><p>下图也展示了差异蛋白富集显著性top 20的GO term，只是横坐标为富集因子（rich factor）。
        Rich factor=（a/b）/（c/d），a表示注释到该term的差异蛋白数，b表示注释到term的差异蛋白总数，c表示该term的背景蛋白数，
        d表示总term的背景蛋白数</p>
        <p>以全部鉴定的蛋白质为背景参照结果目录：<span class="text-primary">6. DEPFunction\*.vs.*.all_GO_Ident\[BP|CC|MF].updown.RFBar</span></p>
        <p>以该物种全部基因为背景参照结果目录:<span class="text-primary">6. DEPFunction\*.vs.*.all_GO_Species\[BP|CC|MF].updown.RFBar</span></p>
        </div>
        """)
        self.print_accordion_head(dpo, "updownrfbar")
        c = 0
        for comp in self.comp_group:
            if c == 0:
                show = 'show'
            else:
                show = ''
            content = self.return_image(
                depfuncurl + '/' + comp + '/' + comp + '.all_GO_Ident/BP.updown.RFBar.png',
                "图{}-{} {}组的差异蛋白的BP term功能富集bar图（rich factor）".format(self.index, figid, comp))
            figid = figid + 1

            content = content + self.return_image(
                depfuncurl + '/' + comp + '/' + comp + '.all_GO_Ident/CC.updown.RFBar.png',
                "图{}-{} {}组的差异蛋白的CC term功能富集bar图（rich factor）".format(self.index, figid, comp))
            figid = figid + 1

            content = content + self.return_image(
                depfuncurl + '/' + comp + '/' + comp + '.all_GO_Ident/MF.updown.RFBar.png',
                "图{}-{} {}组的差异蛋白的MF term功能富集bar图（rich factor）".format(self.index, figid, comp))
            figid = figid + 1

            self.print_accordion_element(dpo, "updownrfbar", "updownrfbar" + str(c),
                                         "updownrf" + str(c), comp,
                                         content, show
                                         )
            figid = figid + 1
            c = c + 1

        self.print_accordion_bottom(dpo)


        #top20 气泡图
        dpo.write("""<div><br /><br /><p>下图展示了差异蛋白富集显著性top 10的GO term，气泡图可以同时显示p-value、rich factor和count。
        横坐标表示rich factor，纵坐标表示p-value的负对数转化，圆圈大小表示count，不同的颜色表示各个term。</p>
                <p>以全部鉴定的蛋白质为背景参照结果目录：<span class="text-primary">6. DEPFunction\*.vs.*.all_GO_Ident\Top10[BP|CC|MF].EnrichedSymbol</span></p>
                <p>以该物种全部基因为背景参照结果目录:<span class="text-primary">6. DEPFunction\*.vs.*.all_GO_Species\[Top10[BP|CC|MF].EnrichedSymbol</span></p>
                </div>
                """)
        self.print_accordion_head(dpo, "top10enrichedsymbol")
        c = 0
        for comp in self.comp_group:
            if c == 0:
                show = 'show'
            else:
                show = ''
            content = self.return_image(
                depfuncurl + '/' + comp + '/' + comp + '.all_GO_Ident/Top10BP.EnrichedSymbol.png',
                "图{}-{} {}组的差异蛋白的BP term富集气泡图（top10）".format(self.index, figid, comp))
            figid = figid + 1

            content = content + self.return_image(
                depfuncurl + '/' + comp + '/' + comp + '.all_GO_Ident/Top10CC.EnrichedSymbol.png',
                "图{}-{} {}组的差异蛋白的CC term富集气泡图（top10）".format(self.index, figid, comp))
            figid = figid + 1

            content = content + self.return_image(
                depfuncurl + '/' + comp + '/' + comp + '.all_GO_Ident/Top10MF.EnrichedSymbol.png',
                "图{}-{} {}组的差异蛋白的MF term富集气泡图（top10）".format(self.index, figid, comp))
            figid = figid + 1

            self.print_accordion_element(dpo, "top10enrichedsymbol", "top10enrichedsymbol" + str(c),
                                         "top10enrsym" + str(c), comp,
                                         content, show
                                         )
            figid = figid + 1
            c = c + 1

        self.print_accordion_bottom(dpo)



        # top20气泡图
        dpo.write("""<div><br /><br /><p>下图展示了差异蛋白富集显著性top20的term。横坐标表示p-value的负对数转化，圆圈颜色表示rich factor，圆圈大小表示count。</p>
                                <p>以全部鉴定的蛋白质为背景参照结果目录：<span class="text-primary">6. DEPFunction\*.vs.*.all_GO_Ident\Top20.[BP|MF|CC].EnrichedSymbol2</span></p>
                                <p>以该物种全部基因为背景参照结果目录:<span class="text-primary">6. DEPFunction\*.vs.*.all_GO_Species\Top20.[BP|MF|CC].EnrichedSymbol2</span></p>
                                </div>
                                """)
        self.print_accordion_head(dpo, "Top20enrichedsymbol")
        c = 0
        for comp in self.comp_group:
            if c == 0:
                show = 'show'
            else:
                show = ''
            content = self.return_image(
                depfuncurl + '/' + comp + '/' + comp + '.all_GO_Ident/Top20.BP.EnrichedSymbol2.png',
                "图{}-{} {}组的差异蛋白的BP term富集气泡图（top20）".format(self.index, figid, comp))
            figid = figid + 1

            content = content + self.return_image(
                depfuncurl + '/' + comp + '/' + comp + '.all_GO_Ident/Top20.CC.EnrichedSymbol2.png',
                "图{}-{} {}组的差异蛋白的CC term富集气泡图（top20）".format(self.index, figid, comp))
            figid = figid + 1

            content = content + self.return_image(
                depfuncurl + '/' + comp + '/' + comp + '.all_GO_Ident/Top20.MF.EnrichedSymbol2.png',
                "图{}-{} {}组的差异蛋白的MF term富集气泡图（top20）".format(self.index, figid, comp))
            figid = figid + 1

            self.print_accordion_element(dpo, "Top20enrichedsymbol", "Top20enrichedsymbol" + str(c),
                                         "top20enrsym" + str(c), comp,
                                         content, show
                                         )
            figid = figid + 1
            c = c + 1

        self.print_accordion_bottom(dpo)

        # DAG图
        dpo.write("""<div><br /><br /><p>下图展示了差异蛋白富集term的有向无循环图。对显著性最高的前10个term作为有向无循环图的主节点，
                以方框表示。并通过包含关系（is a）将相关联的term一起展示，用圆框表示。包含关系（is a）的具体解释可参见
                http://geneontology.org/docs/ontology-relations/。方框或圆框中展示了term的ID、description、p值、count与背景数。
                颜色越接近红色富集显著性越高。</p>
                <p>以全部鉴定的蛋白质为背景参照结果目录：<span class="text-primary">6. DEPFunction\*.vs.*.all_GO_Ident\[BP|MF|CC].DAG.png</span></p>
                <p>以该物种全部基因为背景参照结果目录:<span class="text-primary">6. DEPFunction\*.vs.*.all_GO_Species\[BP|MF|CC].DAG.png</span></p>
                </div>
                                """)
        self.print_accordion_head(dpo, "dagraph")
        c = 0
        for comp in self.comp_group:
            if c == 0:
                show = 'show'
            else:
                show = ''
            content = self.return_image(
                depfuncurl + '/' + comp + '/' + comp + '.all_GO_Ident/BP.DAG.png',
                "图{}-{} {}组的差异蛋白的BP富集有向无循环图".format(self.index, figid, comp))
            figid = figid + 1

            content = content + self.return_image(
                depfuncurl + '/' + comp + '/' + comp + '.all_GO_Ident/CC.DAG.png',
                "图{}-{} {}组的差异蛋白的CC富集有向无循环图".format(self.index, figid, comp))
            figid = figid + 1

            content = content + self.return_image(
                depfuncurl + '/' + comp + '/' + comp + '.all_GO_Ident/MF.DAG.png',
                "图{}-{} {}组的差异蛋白的MF富集有向无循环图".format(self.index, figid, comp))
            figid = figid + 1

            self.print_accordion_element(dpo, "dagraph", "dagraph" + str(c),
                                         "dag" + str(c), comp,
                                         content, show
                                         )
            figid = figid + 1
            c = c + 1

        self.print_accordion_bottom(dpo)

        #KEGG通路分析
        #dpo.write("<h2>{}.{} 差异蛋白的KEGG pathway注释与富集</h2>".format(self.index, h2id))
        self.print_h2(dpo, "差异蛋白的KEGG pathway注释与富集", h2id, "dep.html",
                      "dep{}".format(h2id))
        h2id = h2id + 1

        dpo.write("""<div><p>KEGG数据库（Kyoto Encyclopedia of Genes and Genomes）用于从分子水平的信息中了解细胞、生物体和生态系
        统的高级功能和用途，尤其是从基因组测序和其他高通量实验技术生成的大规模分子数据集中获取信息。在生物体中，蛋白质并不独立行使其功能，而
        是不同蛋白质相互协调完成一系列生化反应以行使其生物学功能。因此，KEGG通路（KEGG pathway）分析是更系统、全面地了解细胞的生物学过程、
        性状或疾病的发生机理、药物作用机制等最直接和必要的途径，将有助于我们去了解蛋白的功能和蛋白间的关系。</p><p>对蛋白进行通路注释，下图
        展示了某个pathway的注释情况，图中的方框为蛋白质/基因，圆点为代谢物，其中红色填充方框表示上调表达的蛋白，黄色填充方框表示下调表达的
        蛋白，绿色方框表示背景蛋白。</p>
        <p>结果目录：6. DEPFunction\\*.vs.*\\*.vs.*\\MX.vs.XZ.all_KEGG\\KEGG.enriched.xlsx，可以通过列表中的URL查看并下载相关通路注释图谱。
</p></div>""")
        dpo.write("""<div><p>依据Fisher精确检验的算法对比较组差异蛋白进行功能富集分析。结果用p-value表示，小于0.05为功能显著富集，
        值越小表示功能富集越显著。下图展示了差异蛋白富集显著性top20的KEGG pathway，横坐标分别为p-value的负对数转化，纵坐标为pathway。
        图中还展示了同一pathway中的上调&下调蛋白数量比例，用不同颜色表示，每个柱子右侧有标注具体的count（富集到该pathway的差异蛋白数量）
        和p值。</p><p>结果目录：<span class="text-primary">6. DEPFunction\\*.vs.*.all_KEGG\\KEGG.enriched.updown.pvalueBar</span></p>
        </div>""")

        self.print_accordion_head(dpo, "keggupdownbar")
        c = 0
        for comp in self.comp_group:
            if c == 0:
                show = 'show'
            else:
                show = ''
            self.print_accordion_element(dpo, "keggupdownbar", "keggupdownbar" + str(c),
                                         "keggupdown" + str(c), comp,
                                         self.return_image(
                                             depfuncurl + '/' + comp + '/' + comp + '.all_KEGG/KEGG.enriched.updown.pvalueBar.png',
                                             "图{}-{} {}组的差异蛋白的KEGG pathway富集bar图（p-value，up&down）".format(self.index, figid, comp)), show
                                         )
            figid = figid + 1
            c = c + 1

        self.print_accordion_bottom(dpo)

        #KEGG RF
        dpo.write("""<div><br /><br /><p>下图也展示了差异蛋白富集显著性top 20的pathway，只是横坐标为富集因子（rich factor）。
        Rich factor=（a/b）/（c/d），a表示注释到该pathway的差异蛋白数，b表示注释到pathway的差异蛋白总数，c表示该pathway的背景蛋白数，
        d表示总pathway的背景蛋白数</p><p>结果目录：<span class="text-primary">
        6. DEPFunction\\*.vs.*.all_KEGG\\KEGG.enriched.updown.RFBar</span></p>
                </div>""")

        self.print_accordion_head(dpo, "keggupdownrfbar")
        c = 0
        for comp in self.comp_group:
            if c == 0:
                show = 'show'
            else:
                show = ''
            self.print_accordion_element(dpo, "keggupdownrfbar", "keggupdownrfbar" + str(c),
                                         "keggupdownrf" + str(c), comp,
                                         self.return_image(
                                             depfuncurl + '/' + comp + '/' + comp + '.all_KEGG/KEGG.enriched.updown.RFBar.png',
                                             "图{}-{} {}组的差异蛋白的KEGG pathway富集bar图（Rich factor，up&down）".format(self.index,
                                                                                                          figid, comp)),
                                         show
                                         )
            figid = figid + 1
            c = c + 1

        self.print_accordion_bottom(dpo)

        # KEGG 气泡图
        dpo.write("""<div><br /><br /><p>下图展示了富集显著性top 10的pathway，用气泡图展示，可以将p-value、Rich factor和count3
        个参数同时体现在一张图内。横坐标表示rich factor，纵坐标表示p-value的负对数转化，圆圈大小表示count，不同的颜色表示各个pathway。
        </p><p>结果目录：<span class="text-primary">
                6. DEPFunction\\*.vs.*.all_KEGG\\Top10EnrichSymbol</span></p>
                        </div>""")

        self.print_accordion_head(dpo, "top10enrichsymbol")
        c = 0
        for comp in self.comp_group:
            if c == 0:
                show = 'show'
            else:
                show = ''
            self.print_accordion_element(dpo, "top10enrichsymbol", "top10enrichsymbol" + str(c),
                                         "top10enrsym" + str(c), comp,
                                         self.return_image(
                                             depfuncurl + '/' + comp + '/' + comp + '.all_KEGG/Top10EnrichSymbol.png',
                                             "图{}-{} {}组的差异蛋白的KEGG通路富集气泡图（top10）".format(
                                                 self.index,
                                                 figid, comp)),
                                         show
                                         )
            figid = figid + 1
            c = c + 1

        self.print_accordion_bottom(dpo)

        # KEGG 气泡图2
        dpo.write("""<div><br /><br /><p>下图展示了显著性top20的pathway，同样用气泡图展示。横坐标表示p-value的负对数转化，圆圈颜色
        表示rich factor，圆圈大小表示count。
                </p><p>结果目录：<span class="text-primary">
                        6. DEPFunction\\*.vs.*.all_KEGG\\Top20EnrichedSymbol2</span></p>
                                </div>""")

        self.print_accordion_head(dpo, "top20enrichsymbol")
        c = 0
        for comp in self.comp_group:
            if c == 0:
                show = 'show'
            else:
                show = ''
            self.print_accordion_element(dpo, "top20enrichsymbol", "top20enrichsymbol" + str(c),
                                         "top20enrbar" + str(c), comp,
                                         self.return_image(
                                             depfuncurl + '/' + comp + '/' + comp + '.all_KEGG/Top20EnrichedSymbol2.png',
                                             "图{}-{} {}组的差异蛋白的KEGG通路富集气泡图（top20）".format(
                                                 self.index,
                                                 figid, comp)),
                                         show
                                         )
            figid = figid + 1
            c = c + 1

        self.print_accordion_bottom(dpo)

        # KEGG 条形图top class
        dpo.write("""<div><br /><br /><p>下图展示了差异蛋白通路富集显著性top30的KEGG pathway。横坐标为p-value的负对数转化，纵坐标
        为具体通路。柱子颜色归属level 1分类。柱子右侧有具体的count和p-value。
                        </p><p>结果目录：<span class="text-primary">
                                6. DEPFunction\\*.vs.*.all_KEGG\\KEGG.Sig.Bar</span></p>
                                        </div>""")

        self.print_accordion_head(dpo, "keggsigbar")
        c = 0
        for comp in self.comp_group:
            if c == 0:
                show = 'show'
            else:
                show = ''
            self.print_accordion_element(dpo, "keggsigbar", "keggsigbar" + str(c),
                                         "keggsig" + str(c), comp,
                                         self.return_image(
                                             depfuncurl + '/' + comp + '/' + comp + '.all_KEGG/KEGG.Sig.Bar.png',
                                             "图{}-{} {}组的差异蛋白的KEGG pathway富集bar图".format(
                                                 self.index,
                                                 figid, comp)),
                                         show
                                         )
            figid = figid + 1
            c = c + 1

        self.print_accordion_bottom(dpo)

        # KEGG 条形图secondary class
        dpo.write("""<div><br /><br /><p>下图也展示了差异蛋白通路富集显著性top30的KEGG pathway。柱子颜色归属level 2分类。
        </p><p>结果目录：<span class="text-primary">
                                6. DEPFunction\\*.vs.*.all_KEGG\\KEGG.Sig.Bar2</span></p>
                                        </div>""")

        self.print_accordion_head(dpo, "keggsig2bar")
        c = 0
        for comp in self.comp_group:
            if c == 0:
                show = 'show'
            else:
                show = ''
            self.print_accordion_element(dpo, "keggsig2bar", "keggsig2bar" + str(c),
                                         "keggsig2" + str(c), comp,
                                         self.return_image(
                                             depfuncurl + '/' + comp + '/' + comp + '.all_KEGG/KEGG.Sig.Bar2.png',
                                             """图{}-{} {}组的差异蛋白的KEGG pathway富集bar图（p-value，level 2 pathway）""".format(
                                                 self.index,
                                                 figid, comp)),
                                         show
                                         )
            figid = figid + 1
            c = c + 1

        self.print_accordion_bottom(dpo)

        # KEGG top class气泡图
        dpo.write("""<div><br /><br /><p>下图也展示了差异蛋白通路富集显著性top 30的pathway，用气泡图展示。将pathway归属到level 1
        分类，每个类的pathway从上至下-log10（p-value）依次降低，即p-value依次升高，显著性依次降低。圆圈大小表示count，颜色表示
        rich factor。
                </p><p>结果目录：<span class="text-primary">
                                        6. DEPFunction\\*.vs.*.all_KEGG\\KEGG.bubble</span></p>
                                                </div>""")

        self.print_accordion_head(dpo, "keggbubble")
        c = 0
        for comp in self.comp_group:
            if c == 0:
                show = 'show'
            else:
                show = ''
            self.print_accordion_element(dpo, "keggbubble", "keggbubble" + str(c),
                                         "keggbub" + str(c), comp,
                                         self.return_image(
                                             depfuncurl + '/' + comp + '/' + comp + '.all_KEGG/KEGG.bubble.png',
                                             """图{}-{} {}组的差异蛋白的KEGG pathway富集气泡图""".format(
                                                 self.index,
                                                 figid, comp)),
                                         show
                                         )
            figid = figid + 1
            c = c + 1

        self.print_accordion_bottom(dpo)

        # KEGG 上下调图
        dpo.write("""<div><br /><br /><p>下图展示了上调蛋白和下调蛋白的通路富集情况（仅all diff的结果中有），横坐标为p-value的对数（下调基因）和负对数（上调基因）转化，X轴大于0，表示上调基因富集的通路，X轴越大p-value越显著；
X轴小于0，表示下调基因富集的通路，X轴越小p-value越显著。蓝色柱子表示下调蛋白富集到的通路，红色柱子表示上调蛋白富集到的通路，分别展示了富集显著性top 12的通路。
                        </p><p>结果目录：<span class="text-primary">
                                                6. DEPFunction\\*.vs.*.all_KEGG\\updownSbar</span></p>
                                                        </div>""")

        self.print_accordion_head(dpo, "keggupdownSbar")
        c = 0
        for comp in self.comp_group:
            if c == 0:
                show = 'show'
            else:
                show = ''
            self.print_accordion_element(dpo, "keggupdownSbar", "keggupdownSbar" + str(c),
                                         "keggupdow" + str(c), comp,
                                         self.return_image(
                                             depfuncurl + '/' + comp + '/' + comp + '.all_KEGG/updownSbar.png',
                                             """图{}-{} {}组的上调&下调蛋白的KEGG通路富集蝴蝶图""".format(
                                                 self.index,
                                                 figid, comp)),
                                         show
                                         )
            figid = figid + 1
            c = c + 1

        self.print_accordion_bottom(dpo)

        #细胞定位分析
        #dpo.write("<h2>{}.{} 差异蛋白的亚细胞定位分析（CC）</h2>".format(self.index, h2id))
        self.print_h2(dpo, "差异蛋白的亚细胞定位分析", h2id, "dep.html",
                      "dep{}".format(h2id))
        h2id = h2id + 1
        dpo.write("""<div><p>细胞可以分为多个细胞器或者细胞区域，如细胞核（nucleolus），细胞质（cytoplasm），线粒体（mitochondrion），
        叶绿体（chloroplast），胞外区（extracellular region），核膜（nuclear membrane），过氧化物酶体（peroxisome），内质网
        （endoplasmic reticulum），高尔基体（Golgi apparatus），溶酶体（lysosome），液泡（vacuole）等，这些细胞器被称为亚细胞。
        不同的基因产物只有在特定的亚细胞中才能正常发挥功能，因此预测蛋白亚细胞定位对了解蛋白结构和特性以及蛋白间的相互作用具有重要意义。
        通过分析GO数据库的细胞组分（Cellular Component，CC），对差异蛋白的亚细胞定位进行注释和统计，下图展示了亚细胞定位的差异蛋白数量
        统计结果，依次为柱形图、气泡图、饼图和环图。下图以全部差异蛋白的亚细胞定位分析结果为例</p>
        <p>结果目录：<span class="text-primary">6. DEPFunction\\*.vs.*\\*.vs.*.all/up/down_subcellular_localization</span></p></div>""")
        self.print_accordion_head(dpo, "subcellularlocalization")
        c = 0
        for comp in self.comp_group:
            if c == 0:
                show = 'show'
            else:
                show = ''
            self.print_accordion_element(dpo, "subcellularlocalization", "subcellularlocalization" + str(c),
                                         "subloc" + str(c), comp,
                                         self.return_image(
                                             depfuncurl + '/' + comp + '/' + comp + '.all_subcellular_localization/CellLocationMerge.png',
                                             """图{}-{} {}组差异蛋白的CC亚细胞定位结果""".format(
                                                 self.index,
                                                 figid, comp)),
                                         show
                                         )
            figid = figid + 1
            c = c + 1

        self.print_accordion_bottom(dpo)

        ##########网络分析###########
        #dpo.write("<h2>{}.{} 蛋白功能互作网络（PFIN）分析</h2>".format(self.index, h2id))
        self.print_h2(dpo, "蛋白功能互作网络（PFIN）分析", h2id, "dep.html",
                      "dep{}".format(h2id))
        h2id = h2id + 1

        dpo.write("""<div><p>在生物体中，蛋白质并不是独立存在的，其功能的行使必须借助于其他蛋白质的调节和介导。这种调节或介导作用的实现
        首先要求蛋白质之间有结合作用或相互作用。对蛋白质之间的相互作用及功能相互作用形成的网络进行研究，对于揭示蛋白质的功能具有重要意义。
        例如，高度聚集的蛋白质可能具有相同或相似的功能；连接度高的蛋白质可能是影响整个系统代谢或信号转导途径的关键点。我们结合来自于STRING数
        据库的蛋白质-蛋白质相互作用关系，以及通路与蛋白之间的关系，整合成一个网络</p>
        <p>结果目录：<span class="text-primary">6.DEPFunction\\*.vs.*\\*.vs.*.PFIN</p>
        </div>""")
        dpo.write("<div>结果链接表：")
        dpo.write("""<table class="table">
            <thead>
              <tr>
                <th>比较组</th>
                <th>网络类型</th>
                <th>结果链接</th>
              </tr>
            </thead>
            <tbody>
            """)
        for comp in self.comp_group:
            if os.path.exists(depfun + '/' + comp + '/' + comp + ".PFIN/pathway-gene_all/web_session/index.html"):
                dpo.write("""
          <tr>
            <td>{0}</td>
            <td>全部显著富集通路(p&lt;0.05)与该通路中蛋白/蛋白编码基因的网络</td>
            <td><a href="{1}" target=_blank><span class="badge-success">结果链接</span></a></td>
          </tr>
          """.format(comp,depfuncurl + '/' + comp + '/' + comp + ".PFIN/pathway-gene_all/web_session/index.html"))
            if os.path.exists(depfun + '/' + comp + '/' + comp + ".PFIN/pathway-gene_part/web_session/index.html"):
                dpo.write("""
                <tr>
                  <td>{0}</td>
                  <td>Top5显著富集通路(p&lt;0.05)与该通路中蛋白/蛋白编码基因的网络</td>
                  <td><a href="{1}" target=_blank><span class="badge-success">结果链接</span></a></td>
                </tr>
                """.format(comp, depfuncurl + '/' + comp + '/' + comp + ".PFIN/pathway-gene_part/web_session/index.html"))
            if os.path.exists(depfun + '/' + comp + '/' + comp + ".PFIN/pathway-gene-ppi_all/web_session/index.html"):
                dpo.write("""
                <tr>
                  <td>{0}</td>
                  <td>Top5显著富集通路(p&lt;0.05)与显著差异蛋白相互作用的网络</td>
                  <td><a href="{1}" target=_blank><span class="badge-success">结果链接</span></a></td>
                </tr>
                """.format(comp, depfuncurl + '/' + comp + '/' + comp + ".PFIN/pathway-gene-ppi_all/web_session/index.html"))
            if os.path.exists(depfun +  '/' + comp + '/' + comp + ".PFIN/pathway-gene-ppi_part/web_session/index.html"):
                dpo.write("""
                <tr>
                  <td>{0}</td>
                  <td>Top5显著富集通路(p&lt;0.05)与该通路中的蛋白/蛋白编码基因及蛋白相互作用网络</td>
                  <td><a href="{1}" target=_blank><span class="badge-success">结果链接</span></a></td>
                </tr>
                """.format(comp, depfuncurl + '/' + comp + '/' + comp + ".PFIN/pathway-gene-ppi_part/web_session/index.html"))
            if os.path.exists(depfun + '/' + comp + '/' + comp + ".PFIN/pathway-pathway/web_session/index.html"):
                dpo.write("""
                   <tr>
                     <td>{0}</td>
                     <td>显著富集通路之间的关联网络(通过共有蛋白关联)</td>
                     <td><a href="{1}" target=_blank><span class="badge-success">结果链接</span></a></td>
                   </tr>
                   """.format(comp,
                              depfuncurl + '/' + comp + '/' + comp + ".PFIN/pathway-pathway/web_session/index.html"))
            if os.path.exists(depfun + '/' + comp + '/' + comp + ".PFIN/ppi/web_session/index.html"):
                dpo.write("""
                              <tr>
                                <td>{0}</td>
                                <td>蛋白相互作用网络</td>
                                <td><a href="{1}" target=_blank><span class="badge-success">结果链接</span></a></td>
                              </tr>
                              """.format(comp,
                                         depfuncurl + '/' + comp + '/' + comp + ".PFIN/ppi/web_session/index.html"))

        dpo.write("""
    </tbody>
</table>""")
        dpo.write("</div><br /><br />")

        self.print_bottom(dpo)
        self.index = self.index + 1


    def print_qc_analysis(self):
        #self.index=index
        qco = open(self.htmlout + '/qc.html', 'w')
        self.print_head(qco)
        qcfold = self.pn + '/2.QualityControl'
        qcurl = '../2.QualityControl'
        figid = 1
        if os.path.exists(qcfold):
            self.print_h1(qco, "数据质控", 'qc.html', "qc")
            qco.write("""<div><p>运用统计学方法，可以直观地反映出搜库结果和数据质量。质控目标有以下三个方面：1. 全，鉴定数目多，覆盖
            率高（主观判断）；2. 稳，鉴定重现性高，稳定性好，尤其是大样本的检测（客观判断）；3. 准，定量准确性高（客观判断）。本实验采用
            高质量精度、高分辨率质谱仪，在数据采集过程中可以保持良好的质量偏差，肽段质量误差的阈值设为±10ppm，结合搜库定性的阈值
            Peptide FDR≤0.01 和 Protein FDR≤0.01 作为筛选标准，以及肽段得分分布等，进一步说明后续用于分析的数据质量较高，说明
            鉴定结果准确可靠。</p></div>
            <p>结果目录：<a href="{}" target=_blank>2.QualityControl</p>
            """.format(qcurl))

        """
        if os.path.exists(qcfold+'/peptide_error.png'):
            self.print_image(qco, qcurl+'/peptide_error.png', "图{}-{}. 肽段误差分布散点图".format(self.index, figid))
            figid = figid + 1
        """

        if os.path.exists(qcfold + '/Peptide_length.png'):
            self.print_image(qco, qcurl+'/Peptide_length.png', "图{}-{}. 肽段长度分布图".format(self.index, figid))
            figid = figid + 1

        if os.path.exists(qcfold + '/unique_peptide.png'):
            self.print_image(qco, qcurl+'/unique_peptide.png', "图{}-{}. unique肽段数量分布图".format(self.index, figid))
            figid = figid + 1

        if os.path.exists(qcfold + '/protein_MW&CI2.png'):
            self.print_image(qco, qcurl+'/protein_MW&CI2.png', "图{}-{}. 蛋白分子量与等电点分布图".format(self.index, figid))
            figid = figid + 1

        if os.path.exists(qcfold + '/protein_sequence_coverage.png'):
            self.print_image(qco, qcurl+'/protein_sequence_coverage.png', "图{}-{}. 蛋白覆盖度分布图".format(self.index, figid))
            figid = figid + 1

        self.print_bottom(qco)
        self.index = self.index + 1

    def print_sample_compare_analysis(self):
        #self.index = index

        scfold = self.pn + '/3.SampleAnalysis'
        scurl = '../3.SampleAnalysis'
        figid = 1

        if os.path.exists(scfold):
            sco = open(self.htmlout + '/sample.html', 'w')
            self.print_head(sco)
            self.print_h1(sco, "样本比较分析","sample.html","sample")
            sco.write("""<div><p>根据蛋白在不同样本间的表达情况，对样本进行表达分布、相关性分析和PCA分析。根据结果可以从各个角度评估
            组间差异性和组内重复性。</p><p>以下均为样本表达分布图。箱线图和小提琴图的同组样本中位数接近于同一水平线，直方图接近于钟形分
            布，说明数据质量越好，重复性良好。</p><p>结果目录：<a href="{}" target=_blank>3.SampleAnalysis</a></p>
            </div>""".format(scurl))

            if os.path.exists(scfold+'/boxplot.png'):
                self.print_image(sco, scurl+'/boxplot.png', "图{}-{}. 样本表达分布图的箱线图".format(self.index, figid))

            if os.path.exists(scfold+'/violin.png'):
                self.print_image(sco, scurl+'/violin.png', "图{}-{}. 样本表达分布图的小提琴图".format(self.index, figid))

            if os.path.exists(scfold+'/histogramFacet.png'):
                self.print_image(sco, scurl+'/histogramFacet.png', "图{}-{}. 样本表达分布图的直方图".format(self.index, figid))

            if os.path.exists(scfold+'/corr.png'):
                sco.write("""以下为样本相关性图。相关性系数由-1至1，＜0为负相关，＞0为正相关，绝对值越大表示相关性越强。也就是同组样本间的
                相关性系数越接近于1越好。""")
                self.print_image(sco, scurl+'/corr.png', "图{}-{}. 样本相关性图（气泡图）".format(self.index, figid))

            if os.path.exists(scfold + '/ggpairs.png'):
                self.print_image(sco, scurl + '/ggpairs.png', "图{}-{}. 样本相关性图（散点图）".format(self.index, figid))

            if os.path.exists(scfold + '/PCA1.2.png'):
                sco.write("""以下为样本PCA得分图。每个点表示一个样本，样本间距离越近表明相似性越高，距离越远表示差异性越大。
                背景椭圆表示置信区间，超出置信圈的样本可能为异常样本。""")
                self.print_image(sco, scurl + '/PCA1.2.png', "图{}-{}. 样本PCA得分图".format(self.index, figid))

            if os.path.exists(scfold):
                self.print_bottom(sco)

        self.index = self.index + 1


    def print_multi_group_comparison_analysis(self, index):
        self.index = index
        anovafold = self.pn+'/7.ANOVA'
        anovaurl = '../4.AllProtein/4.1Multi_group_expression/kmeans'
        if os.path.exists(anovafold):
            anovagroup = os.listdir(anovafold)
            mco = open(self.htmlout + '/anova.html', 'w')
            self.print_head(mco)
            self.print_h1(mco, "多组表达分析", 'anova.html', 'multicomp')
            h2id = 1
            figid = 1

            self.print_h2(mco,"单因素方差检验和比较绘图", h2id, 'anova.html', "multicomp{}".format(h2id))
            h2id = h2id + 1

            mco.write("""<div><p>分析蛋白在多组中的表达量变化趋势，采用单因素方差分析（ANOVA）检验方法。<p></div>""")
            if os.path.exists(self.pn + '/7.ANOVA/anova_sig_number.txt'):
                mco.write("""<div><p>单因素方检验中统计学上显著的蛋白质数的统计(p&lt;0.05)</p></div>""")
                mco.write("""<table class="table"><thead><tr><td>组别</td><td>显著的蛋白数</td></tr></thead><tbody>""")
                ias = open(self.pn + '/7.ANOVA/anova_sig_number.txt')
                for line in ias.readlines():
                    row = line.rstrip().split('\t')
                    mco.write("""<tr><td>{}</td><td>{}</td></tr>""".format(row[0], row[1]))
                ias.close()
                mco.write("""</tbody></table>""")


            mco.write("""<div><p>对显著性p-value&lt;0.05的蛋白进行可视化。
            横坐标为不同组，纵坐标为相对表达量。其中BarComp的两组间显著性采用T检验，用*表示，*表示p&lt;0.05，**表示p&lt;0.01，
            ***表示p&lt;0.001，ns表示p&gt;0.05（即没有显著性）。</p>
            <p>注：仅展示显著性top300的蛋白</p><p>
            <p>ANOVA检验结果：<span class="text-primary">
            7.ANOVA\\[group name]\\kmeans\\ANOVA_result.xlsx</span></p></div>
            """)

            mco.write("""<table class="table"><thead><tr><td>组别</td><td>类型</td><td>目录</td></tr></thead>
            <tbody>
            """)
            for at in anovagroup:
                if os.path.exists(anovafold+ '/'+at+'/Multi_group_expression'):
                    if os.path.exists(anovafold+ '/'+at+'/Multi_group_expression/BarComp'):
                        mco.write("""<tr><td>{0}</td><td>{1}</td><td><a href="{2}" 
                        target="_blank">{3}</a></td></tr>""".format(at, '比较柱状图', '../7.ANOVA/'+at+'/Multi_group_expression/BarComp',"BarComp"))
                    if os.path.exists(anovafold + '/' + at + '/Multi_group_expression/BarSE'):
                        mco.write("""<tr><td>{0}</td><td>{1}</td><td><a href="{2}" 
                                                target="_blank">{3}</a></td></tr>""".format(at, '柱状图',
                                                                                            '../7.ANOVA/' + at + '/Multi_group_expression/BarSE',
                                                                                            "BarSE"))
                    if os.path.exists(anovafold + '/' + at + '/Multi_group_expression/Boxplot'):
                        mco.write("""<tr><td>{0}</td><td>{1}</td><td><a href="{2}" 
                                                target="_blank">{3}</a></td></tr>""".format(at, '箱线图',
                                                                                            '../7.ANOVA/' + at + '/Multi_group_expression/Boxplot',
                                                                                            "Boxplot"))

                    if os.path.exists(anovafold + '/' + at + '/Multi_group_expression/Jitter'):
                        mco.write("""<tr><td>{0}</td><td>{1}</td><td><a href="{2}" 
                                                target="_blank">{3}</a></td></tr>""".format(at, '抖图',
                                                                                            '../7.ANOVA/' + at + '/Multi_group_expression/Jitter',
                                                                                            "Jitter"))
                    if os.path.exists(anovafold + '/' + at + '/Multi_group_expression/ScatterSE'):
                        mco.write("""<tr><td>{0}</td><td>{1}</td><td><a href="{2}" 
                                                target="_blank">{3}</a></td></tr>""".format(at, '散点图',
                                                                                            '../7.ANOVA/' + at + '/Multi_group_expression/ScatterSE',
                                                                                            "ScatterSE"))
            mco.write("""</tbody></table>""")

            #mco.write("<h2>{} k-means聚类分析</h2>".format(str(self.index) + '.' + str(h2id)))
            self.print_h2(mco, "k-means聚类分析", h2id, 'anova.html',
                          "multicomp{}".format(h2id))
            h2id = h2id + 1
            mco.write("""<div><p>聚类算法是指将一堆没有标签的数据自动划分成几类的方法，这个方法要保证同一类的数据有相似的特征。而k均值
            聚类算法（k-means clustering algorithm）是一种迭代求解的聚类分析算法，其步骤是，预将数据分为K组，则随机选取K个对象作为初
            始的聚类中心，然后计算每个对象与各个种子聚类中心之间的距离，把每个对象分配给距离它最近的聚类中心。聚类中心以及分配给它们的对象
            就代表一个聚类。我们将ANOVA检验结果中显著的蛋白（p&lt;0.05)的表达值做了Z-score变换后然使用Hartigan-Wong方法进行k-means聚类分析。</p>
            <p>下图是全部样本的k-means聚类图，绿色线代表该类中每一个蛋白的表达趋势，红色线表示该类的均值。</p></div>""")
            anovaurl="../7.ANOVA"
            self.print_accordion_head(mco , "kmeans")
            c = 0
            for comp in anovagroup:
                if os.path.exists(anovafold+ '/'+comp+'/Multi_group_expression'):
                    if c == 0:
                        show = 'show'
                    else:
                        show = ''

                    content = self.return_image(anovaurl + '/' + comp + '/kmeans/zscoreHartigan-Wong.kmeans.png',
                                                                   "图{}-{} {}组的k-means聚类分析图".format(self.index, figid, comp))+ self.return_image(anovaurl + '/' + comp + '/kmeans/zscoreavgHartigan-Wong.kmeans.png',
                                      "图{}-{} {}组的每组均值的k-means聚类分析图".format(self.index, figid, comp))

                    self.print_accordion_element(mco, "kmeans", "kmeans" + str(c),
                                                 "vol" + str(c), comp,
                                                 content, show
                                                 )
                    figid = figid + 1
                    c = c + 1

            self.print_accordion_bottom(mco)
            self.print_h2(mco, "GO功能和KEGG通路分析", h2id, 'anova.html',
                          "multicomp{}".format(h2id))
            h2id = h2id + 1
            mco.write("<div><p>ANOVA检验得到的p&lt;0.05的蛋白质的功能分析结果说明可以参考差异蛋白的功能分析结果:</p></div>")
            mco.write("""<table class="table"><thead><tr><td>组别</td><td>类型</td><td>目录</td></tr></thead>
            <tbody>""")

            for at in anovagroup:
                if os.path.exists(anovafold+ '/'+at):
                    if os.path.exists(anovafold+ '/'+at+'/GOEnrich_Identified'):
                        mco.write("""<tr><td>{0}</td><td>{1}</td><td><a href="{2}" 
                        target="_blank">{3}</a></td></tr>""".format(at, '以全部鉴定到蛋白质为背景参照的GO富集分析结果', '../7.ANOVA/'+at+'/GOEnrich_Identified','7.ANOVA/'+at+'/GOEnrich_Identified'))
                    if os.path.exists(anovafold + '/' + at + '/GOEnrich_Species'):
                        mco.write("""<tr><td>{0}</td><td>{1}</td><td><a href="{2}" 
                                                target="_blank">{3}</a></td></tr>""".format(at, '以该物种全部基因为背景参照的GO富集分析结果',
                                                                                            '../7.ANOVA/'+at+'/GOEnrich_Species',
                                                                                            '7.ANOVA/'+at+'/GOEnrich_Species'))
                    if os.path.exists(anovafold + '/' + at + '/GOLevel2'):
                        mco.write("""<tr><td>{0}</td><td>{1}</td><td><a href="{2}" 
                                                target="_blank">{3}</a></td></tr>""".format(at,
                                                                                            'GO第二层功能统计结果',
                                                                                            '../7.ANOVA/' + at + '/GOLevel2',
                                                                                            '7.ANOVA/' + at + '/GOLevel2'))
                    if os.path.exists(anovafold + '/' + at + '/GO_subcellular_localization'):
                        mco.write("""<tr><td>{0}</td><td>{1}</td><td><a href="{2}" 
                                                target="_blank">{3}</a></td></tr>""".format(at,
                                                                                            '蛋白质在经典的亚细胞器定位统计结果',
                                                                                            '../7.ANOVA/' + at + '/GO_subcellular_localization',
                                                                                            '7.ANOVA/' + at + '/GO_subcellular_localization'))
                    if os.path.exists(anovafold + '/' + at + '/KEGG'):
                        mco.write("""<tr><td>{0}</td><td>{1}</td><td><a href="{2}" 
                                                target="_blank">{3}</a></td></tr>""".format(at,
                                                                                            'KEGG通路分析结果',
                                                                                            '../7.ANOVA/' + at + '/KEGG',
                                                                                            '7.ANOVA/' + at + '/KEGG'))
            mco.write("</tbody></table>")

            self.print_bottom(mco)
            self.index = self.index + 1

    def print_advance_html(self):
        ado = open(self.pn+'/html/advance.html', 'w')
        self.print_head(ado)
        ado.write("""
<div id="advance"><h1><div><a name="advance"></a>{0}. 高级分析</div></h1></div>

<h2>{0}.1 WGCNA分析<a name="advance1"></a></h2>

<div><p>加权基因共表达网络分析（Weighted gene co-expression network analysis，WGCNA），是一种应用广泛的数据挖掘方法，在蛋白组
分析中，它可以通过计算蛋白之间的表达相关性，将具有表达相关性的蛋白聚类到一个模块中，然后再分析模块与样本特征（包括临床特征、手术方式、治
疗方法等等）之间的相关性，WGCNA搭建了一座样本特征与蛋白表达变化之间的桥梁。WGCNA的一般流程：1. 利用蛋白之间的相互作用模式构建一个网络；
2. 进行模块识别；3. 将模块与外部信息关联；4. 鉴定保守模块和特异模块；5. 找到网络连接度较高的关键蛋白（hub protein)，从而筛选出目标蛋白。
WGCNA在组学中的应用场景非常广泛，样品数量建议至少5组3重复进行WGCNA分析有生物学意义。</p></div>

<div class="text-center">
      <figure class="figure">
           <a href="img/图8-1. WGCNA分析.png" class="thumbnail" target=_blank>
            <img src="img/图8-1. WGCNA分析.png" class="figure-img img-fluid d-block mx-auto" width=600 alt="无该结果"  />
            </a>
      <figcaption class="figure-caption text-center">图{0}-1. WGCNA分析</figcaption>
        </figure>
</div>
		<h2>{0}.2 机器学习和图形构建<a name="advance2"></a></h2>
		<div><p>早期蛋白质生物标志物的筛选工作中，研究者倾向于获得一种特定的蛋白质作为生物标志物（biomarker），
		但最近有研究表明将多个现有标志物进行组合形成新的评判指标能有效提升预测准确性。此外，将多组学数据综合分析寻找标志物也逐渐成为研究热
		点，多组学来源的生物标志物的组合相较单一来源的生物标志物更为准确。这里可以通过常规的统计方法T-test、FC等对蛋白进行初步过滤，去除一
		些无贡献意义的变量；然后采用ROC分析来筛选出一些候选变量；再采用随机森林、支持向量机等机器学习方法进行建模，通过重要性、表达水平、相
		关性等指标对候选变量进行评价；最后建立诊断panel模型并对其进行评价。本分析建议样品数量至少2组10重复。</p>
                </div><div class="text-center">
                <figure class="figure">
                        <a href="img/图8-2. ROC曲线.png" class="thumbnail" target=_blank>
                            <img src="img/图8-2. ROC曲线.png" class="figure-img img-fluid d-block mx-auto" width=600 
                            alt="无该结果"  />
                        </a>
                        <figcaption class="figure-caption text-center">图{0}-2. ROC曲线（count）</figcaption>
        </figure>
		</div>
		
		<h2>{0}.3 代谢通路整体趋势分析<a name="advance3"></a></h2>
		<div><p>为了系统研究代谢变化，利用差异蛋白丰度对代谢通路进行分析。差异丰度得分（Differential abundance scores，DA score）
		可以捕捉到通路中蛋白相对于对照组增加/减少的趋势。对差异蛋白的注释结果进行DA score的计算，取top 30的通路绘制图片，如下图所示。
		横坐标为DA score，纵坐标为通路，不同颜色归属level 1分类。每个圈对应一条通路，圆圈大小表示注释到该通路中的蛋白数量，颜色由蓝至红
		表示DA score由-1到1。当DA score为-1时，表示该通路中所有蛋白的丰度均下降，当DA score为1时，表示该通路中所有蛋白的丰度均上升。
		</p></div>
                <div class="text-center">
                <figure class="figure">
                        <a href="img/图8-3. 代谢通路整体趋势分析.png" class="thumbnail" target=_blank>
                            <img src="img/图8-3. 代谢通路整体趋势分析.png" class="figure-img img-fluid d-block mx-auto" 
                            width=600 alt="无该结果"  />
                        </a>
                        <figcaption class="figure-caption text-center">图{0}-3. 代谢通路整体趋势分析（rich  factor）</figcaption>
        </figure></div>

		<h2>{0}.4 GSEA蛋白集富集分析(人、其他物种需要定制)<a name="advance4"></a></h2>
		<div><p>基因集富集分析(Gene Set Enrichment Analysis,GSEA)是一种针对全蛋白组表达谱数据的分析方法，将蛋白与预定义的蛋白集进行
		比较。与GO/KEGG富集不同，GO/KEGG富集是先筛选差异蛋白，再判断哪些注释到的GO term / pathway有显著富集，它的结果与差异蛋白的个数
		有关，进一步追溯与差异蛋白的筛选参数的设定有关，存在一定主观性而且只能用于表达变化较大的蛋白，而GSEA不局限于差异蛋白。而且GO/KEGG
		富集筛选到某条通路，这条通路中既有上调蛋白又有下调蛋白，那么这条通路是被抑制还是激活，GSEA可以解释这个问题。</p><p>某通路的GSEA分
		析结果如下图所示，绿色折线表示蛋白的Enrichment score，蛋白从左到右是按照FC由大到小排序，当ES值小于0，该功能富集在蛋白排序的后端，
		表示蛋白低表达，该功能被抑制，反之则激活</p></div>
                <div class="text-center">
                <figure class="figure">
                        <a href="img/图8-4. GSEA分析KEGG pathway(1).png" class="thumbnail" target=_blank>
                            <img src="img/图8-4. GSEA分析KEGG pathway(1).png" class="figure-img img-fluid d-block mx-auto" 
                            width=600 alt="无该结果"  />
                        </a>
						<figure class="figure">
                        <a href="img/图8-4. GSEA分析KEGG pathway(2).png" class="thumbnail" target=_blank>
                            <img src="img/图8-4. GSEA分析KEGG pathway(2).png" class="figure-img img-fluid d-block mx-auto" 
                            width=600 alt="无该结果"  />
                        </a>
                        <figcaption class="figure-caption text-center">图{0}-4. GSEA分析KEGG pathway（rich  factor）</figcaption>
        </figure></div>
		
		<h2>{0}.5 SPIA信号通路富集分析（人、大鼠、小鼠）<a name="advance5"></a></h2>
		<div><p>信号通路影响分析（Signaling Pathway Impact Analysis，SPIA）顾名思义，它是一种分析信号通路的方法。该分析整合了
		FC、p-value和拓扑学分析，为每条通路计算一个总体概率值pG。SPIA考虑了该蛋白在通路中的拓扑学关系，比如差异蛋白之间是上游或下游
		以及激活或抑制的关系。下图富集结果可以显示信号通路激活或抑制的总体概率值pG。GO等注释数据库中不提供拓扑学关系，所以无法对
		GO term进行SPIA富集分析。</p>
                </div><div class="text-center">
                <figure class="figure">
                     <figure class="figure">
                        <a href="img/图8-5. SPIA信号通路富集分析.png" class="thumbnail" target=_blank>
                            <img src="img/图8-5. SPIA信号通路富集分析.png" class="figure-img img-fluid d-block mx-auto" width=600 alt="无该结果"/>
							
                        </a>
						<figcaption class="figure-caption text-center">图{0}-5. SPIA信号通路富集分析</figcaption>
        </figure></div>

		<h2>{0}.6 Wikipathway通路富集分析<a name="advance6"></a></h2>
		<div><p>Wikipathway（https://www.wikipathways.org）是代谢通路专用数据库，支持的物种如下图所示：</p>
                </div>
				<div class="text-center">
                <figure class="figure">
                     <figure class="figure">
                        <a href="img/图8-6. Wikipathway数据库支持的物种.png" class="thumbnail" target=_blank>
                            <img src="img/图8-6. Wikipathway数据库支持的物种.png" class="figure-img img-fluid d-block mx-auto" width=300 alt="无该结果"/>
                        </a>
						<figcaption class="figure-caption text-center">图{0}-6. Wikipathway数据库支持的物种</figcaption>

						<p>Wikipathway结果示例图可以参考上述差异蛋白的KEGG pathway富集，只是两者所用数据库不同。</p>
						<div class="text-center">
                <figure class="figure">
                     <figure class="figure">
                        <a href="img/图8-7. Wikipathway通路富集气泡图.png" class="thumbnail" target=_blank>
                            <img src="img/图8-7. Wikipathway通路富集气泡图.png" class="figure-img img-fluid d-block mx-auto" width=600 alt="无该结果"/>
                        </a>
						<figcaption class="figure-caption text-center">图{0}-7. Wikipathway通路富集气泡图</figcaption>
        </figure></div></div>
		
		<h2>{0}.7 Reactome通路富集分析（人、大鼠、小鼠）<a name="advance7"></a></h2>
		<div><p>REACTOME（https://reactome.org/）数据模型的核心单元是反应。参与反应的实体（核酸，蛋白质，复合物，疫苗，抗癌治疗剂和小分子）形成生物相互作用网络，并组合为通路。Reactome结果示例图可以参考差异蛋白的KEGG pathway富集，只是所用数据库不同。</p>
                </div><div class="text-center">
                <figure class="figure">
                        <a href="img/图8-8. Reactome通路富集分析.png" class="thumbnail" target=_blank>
                            <img src="img/图8-8. Reactome通路富集分析.png" class="figure-img img-fluid d-block mx-auto" width=600 alt="无该结果"  />
                        </a>
                        <figcaption class="figure-caption text-center">图{0}-8. Reactome通路富集分析</figcaption>
        </figure></div>
        """.format(self.index))
        self.print_bottom(ado)
        url = 'advance.html'
        name = 'advance'
        self.lefthandle.write("""<div class="level1"><a href="{}#{}" target="rightframe">{}</div>""".format(url, name,
                                                                                                            str(
                                                                                                                self.index) + '. ' + "高级分析"))
        self.lefthandle.write(
            """<div class="level2"><a href="{}#{}" target="rightframe">{}</div>""".format(url, name + '1',
                                                                                          str(
                                                                                              self.index) + ".1 " + "WGCNA分析"))
        self.lefthandle.write(
            """<div class="level2"><a href="{}#{}" target="rightframe">{}</div>""".format(url, name + '2',
                                                                                          str(
                                                                                              self.index) + ".2 " + "机器学习和图形构建"))
        self.lefthandle.write(
            """<div class="level2"><a href="{}#{}" target="rightframe">{}</div>""".format(url, name + '3',
                                                                                          str(
                                                                                              self.index) + ".3 " + "代谢通路整体趋势分析"))
        self.lefthandle.write(
            """<div class="level2"><a href="{}#{}" target="rightframe">{}</div>""".format(url, name + '4',
                                                                                          str(
                                                                                              self.index) + ".4 " + "GSEA蛋白集富集分析"))
        self.lefthandle.write(
            """<div class="level2"><a href="{}#{}" target="rightframe">{}</div>""".format(url, name + '5',
                                                                                          str(
                                                                                              self.index) + ".5 " + "SPIA信号通路富集分析"))
        self.lefthandle.write(
            """<div class="level2"><a href="{}#{}" target="rightframe">{}</div>""".format(url, name + '6',
                                                                                          str(
                                                                                              self.index) + ".6 " + "Wikipathway通路富集分析"))
        self.lefthandle.write(
            """<div class="level2"><a href="{}#{}" target="rightframe">{}</div>""".format(url, name + '7',
                                                                                          str(
                                                                                              self.index) + ".7 " + "Reactome通路富集分析"))
        ado.close()
        self.index = self.index + 1

    def print_personal_analysis(self):
        po = open(self.pn + '/html/personalized.html', 'w')
        self.print_head(po)
        po.write("""<div id="personal"><h1><div>{0}. 个性化分析(demo展示，非本项目结果，需要请另询！)<a name="personal"></a>
        </div></h1></div>
<h2>{0}.1 复杂热图（complex heatmap）<a name="personal1"></a></h2>
<div><p>常见的简单热图（或称单个热图）只能展示一种元素，比如相对表达量，而复杂热图（complex heatmap）在简单热图的基础上又增加了多个元素，
可以有效地可视化不同数据集之间的关联并揭示其潜在模式。如下图所示，图中展示了相对表达量、FC、p-value和unique肽段数和等多个元素的信息。在数
据分析过程中，老师可以根据需要选择不同的元素进行complex heatmap的定制。</p></div>

<div class="text-center">
      <figure class="figure">
           <a href="img/图9-1. 复杂热图（complex heatmap）图例.png" class="thumbnail" target=_blank>
            <img src="img/图9-1. 复杂热图（complex heatmap）图例.png" class="figure-img img-fluid d-block mx-auto" 
            width=800 alt="无该结果"  />
            </a>
      <figcaption class="figure-caption text-center">图{0}-1. 复杂热图（complex heatmap）图例</figcaption>
        </figure>
</div>
		
		<h2>{0}.2 多组比较数据的富集分析<a name="personal2"></a></h2>
		<div><p>标准分析中提供了比较组的差异代谢物通路富集气泡图，这里可以将各个比较组的通路富集结果整合在一起，以便比较分析。如下图所示，
		横坐标为比较组，纵坐标每个圆圈表示一条通路，圆圈大小表示富集到该通路中的蛋白数量，圆圈颜色表示该通路的富集显著性。</p>
                </div><div class="text-center">
                <figure class="figure">
                        <a href="img图9-2. 多组比较数据通路富集气泡图.png" class="thumbnail" target=_blank>
                            <img src="img/图9-2. 多组比较数据通路富集气泡图.png" class="figure-img img-fluid d-block mx-auto" 
                            width=800 alt="无该结果"  />
                        </a>
                        <figcaption class="figure-caption text-center">图{0}-2. 多组比较数据通路富集气泡图</figcaption>
        </figure>
		</div>
		
		<h2>{0}.3 Circos圈图<a name="personal3"></a></h2>
		<div><p>Circos圈图主要用于数据关系的可视化，其中很重要的一种圈图——弦图在生物学高分文章中经常见到，主要用于展示多个对象之间的关系。
		连接圆上任意两点的弧线叫做弦，表示两点之间存在关联关系，这种形式非常适合分析复杂数据之间的关联关系。另外，热图也可以通过circos圈图
		来展示，相比常规热图，在展示较多的features时圈图的可视化效果更好。</p></div>
                <div class="text-center">
                <figure class="figure">
                        <a href="img/图9-3. Circos图表示蛋白与KEGG pathway之间的关系.png" class="thumbnail" target=_blank>
                            <img src="img/图9-3. Circos图表示蛋白与KEGG pathway之间的关系.png" class="figure-img img-fluid
                             d-block mx-auto" width=800 alt="无该结果"  />
                        </a>
                        <figcaption class="figure-caption text-center">图{0}-3. Circos图表示蛋白与KEGG pathway之间的关系
                        </figcaption>
        </figure></div>
		
				<div class="text-center">
                <figure class="figure">
                        <a href="img/图9-4. Circos图表示蛋白表达热图.png" class="thumbnail" target=_blank>
                            <img src="img/图9-4. Circos图表示蛋白表达热图.png" class="figure-img img-fluid d-block mx-auto"
                             width=800 alt="无该结果"  />
                        </a>
                        <figcaption class="figure-caption text-center">图{0}-3. Circos图表示蛋白与KEGG pathway之间的关系
                        </figcaption>
        </figure></div>

		<h2>{0}.4 桑基图<a name="personal4"></a></h2>
		<div><p>桑基图，又称为桑基能量分流图或桑基能量平衡图，它是一种特定类型的流程图，图中分支曲线的宽度对应数据流量的大小，所有主支宽度
		的总和与所有分支宽度的总和相等，保持能量平衡，适用于可视化分析上下调蛋白与功能或者修饰位点（修饰组学）与蛋白和功能之间的数据流走向
		分析。</p></div>
                <div class="text-center">
                <figure class="figure">
                        <a href="img/图9-5. 桑基图展示上下调蛋白与功能层级的数据流走向分析.png" class="thumbnail" target=_blank>
                            <img src="img/图9-5. 桑基图展示上下调蛋白与功能层级的数据流走向分析.png" class="figure-img
                             img-fluid d-block mx-auto" width=800 alt="无该结果"  />
                        </a>
						
                        <figcaption class="figure-caption text-center">图{0}-5. 桑基图展示上下调蛋白与功能层级的数据流走向分析
                        </figcaption>
        </figure></div>
        """.format(self.index))
        self.print_bottom(po)

        url = 'personalized.html'
        name = 'personal'
        self.lefthandle.write("""<div class="level1"><a href="{}#{}" target="rightframe">{}</div>""".format(url, name,
                                                                                                            str(
                                                                                                                self.index) + '. ' + "个性化分析"))
        self.lefthandle.write(
            """<div class="level2"><a href="{}#{}" target="rightframe">{}</div>""".format(url, name + '1',
                                                                                          str(
                                                                                              self.index) + ".1 " + "复杂热图（complex heatmap）"))
        self.lefthandle.write(
            """<div class="level2"><a href="{}#{}" target="rightframe">{}</div>""".format(url, name + '2',
                                                                                          str(
                                                                                              self.index) + ".2 " + "多组比较数据的富集分析"))
        self.lefthandle.write(
            """<div class="level2"><a href="{}#{}" target="rightframe">{}</div>""".format(url, name + '3',
                                                                                          str(
                                                                                              self.index) + ".3 " + "Circos圈图"))
        self.lefthandle.write(
            """<div class="level2"><a href="{}#{}" target="rightframe">{}</div>""".format(url, name + '4',
                                                                                          str(
                                                                                              self.index) + ".4 " + "桑基图"))

        po.close()
        self.index = self.index + 1

    def print_english_reference(self):
        ero = open(self.pn+'/html/literature.html', 'w')
        self.print_head(ero)
        ero.write("""
	<div id="english">
		<h1>
			<div>{0}. 英文写作参考<a name="lit1"></a></div>
		</h1>
	</div>	
	<div class="alert alert-warning">注：下文仅作为英文模板参考，具体细节要以中文报告为准！</div>
	<h2>Materials and Methods</h2>
	<h2>Sample preparation</h2>
	<p>For cell or tissue lysis, cell pellets /tissue were suspended on ice in 200 μL lysis buffer (4% SDS, 100 mM DTT,
	 150 mM Tris-HCl pH 8.0). Cells /tissue were disrupted with agitation using a homogenizer, and boiling for 5min. 
	 The samples were further ultrasonicated and boiling again for another 5 min. Undissolved cellular debris were 
	 removed by centrifugation at 16000 rpm for 15 min. The supernatant were collected and quantified with a BCA 
	 Protein Assay Kit (Bio-Rad, USA).</p>

	<h2>Protein Digestion</h2>
		<p>Digestion of protein (200 μg for each sample) was performed according to the FASP procedure described by 
		Wisniewski, Zougman et al. Briefly, the detergent, DTT and other low-molecular-weight components were removed
		 using 200 μl UA buffer (8 M Urea, 150 mM Tris-HCl pH 8.0) by repeated ultrafiltration (Microcon units, 30 kD) 
		 facilitated by centrifugation. Then 100 μL 0.05 M iodoacetamide in UA buffer was added to block reduced 
		 cysteine residues and the samples were incubated for 20 min in darkness. The filter was washed with 100 μl UA 
		 buffer three times and then 100 μl 25 mM NH4HCO3 twice. Finally, the protein suspension was digested with 4 μg 
		 trypsin (Promega) in 40 μl 25 mM NH4HCO3 overnight at 37 °C, and the resulting peptides were collected as a 
		 filtrate. The peptide concentration was determined with OD280 by Nanodrop device.</p>
	<h2>TMT Labeling of peptides</h2>
	<p>Peptides were labeled with TMT reagents according to the manufacturer’s instructions (Thermo Fisher Scientific).
	 Each aliquot (100 μg of peptide equivalent) was reacted with one tube of TMT reagent, respectively. After the sample 
	 was dissolved in 100 μL of 0.05M TEAB solution, pH 8.5, the TMT reagent was dissolved in 41 μL of anhydrous 
	 acetonitrile. The mixture was incubated at room temperature for 1 h. Then 8μL of 5% hydroxylamine to the sample and
	  incubate for 15 minutes to quench the reaction. The Multiplex labeled samples were pooled together and 
	  lyophilized.</p>	
	<h2>High pH Reverse Phase Fractionation (HPRP)</h2>
    <p>TMT-labeled peptides mixture was fractionated using a Waters XBridge BEH130 column (C18, 3.5μm, 2.1 × 150mm) on a
     Agilent 1290 HPLC operating at 0.3 mL/min. Buffer A consisted of 10mM ammonium formate and buffer B consisted of 
     10mM ammonium formate with 90% acetonitrile; both buffers were adjusted to pH 10 with ammonium hydroxide. A total
      of 30 fractions were collected for each peptides mixture, and then concatenated to 15 (pooling equal interval 
      RPLC fractions). The fractions were dried for nano LC-MS/MS analysis.</p>	
    <h2>*LC-MS Analysis (TMT10plex)</h2>
    <p>LC- MS analysis were performed on a Q Exactive mass spectrometer that was coupled to Easy nLC (Thermo Fisher
     Scientific). Peptide from each fraction was loaded onto a the C18-reversed phase column (12cm long, 75μm ID, 3μm) 
     in buffer A (2% acetonitrile and 0.1% Formic acid) and separated with a linear gradient of buffer B (90% acetonitrile
      and 0.1% Formic acid) at a flow rate of 300 nL/min over 90 min. The linear gradient was set as follows: 0–2 min, 
      linear gradient from 2% to 5% buffer B; 2–62 min, linear gradient from 5% to 20% buffer B; 62–80 min, linear 
      gradient from 20% to 35% buffer B; 80–83 min, linear gradient from 35% to 90% buffer B; 83–90 min, buffer B maintained
       at 90%. MS data was acquired using a data-dependent top15 method dynamically choosing the most abundant precursor
        ions from the survey scan (300–1800 m/z) for HCD fragmentation. Determination of the target value is based on 
        predictive Automatic Gain Control (pAGC). The AGC target values of 1e6, and maximum injection time 50 ms were 
        for full MS, and a target AGC value of 1e5, maximum injection time 100 ms for MS2. Dynamic exclusion duration 
        was 30s. Survey scans were acquired at a resolution of 70,000 at m/z 200 and resolution for HCD spectra was set
         to 35,000 at m/z 200. Normalized collision energy was 30. The instrument was run with peptide recognition
          mode enabled.</p>
     <h2>Database Searching and Analysis</h2>   
     <p>The resulting LC-MS/MS raw files were imported into Proteome Discoverer 2.4  software (version 1.6.0.16) for 
     data interpretation and protein identification against the database  . (请根据数据库情况进行修改) An initial search
      was set at a precursor mass window of 10 ppm. The search followed an enzymatic cleavage rule of Trypsin/P and 
      allowed maximal two missed cleavage sites and a mass tolerance of 20ppm for fragment ions. The modification set
       was as following: fixed modification: Carbamidomethyl (C), TMT10plex(K), TMT10plex(N-term), Variable 
       modification：Oxidation(M) and Acetyl (Protein N-term). The minimum 6 amino acids for peptide, ≥1 unique peptides
        were required per protein. For peptide and protein identification, false discovery rate (FDR) was set to 1%. TMT
         reporter ion intensity were used for quantification.</p>   
     <h2>Bioinformatics analysis</h2>
     <p>Analyses of bioinformatics data were carried out with Perseus software, Microsoft Excel and R statistical 
     computing software. Differentially significant expressed proteins were screened with the cutoff of a ratio 
     fold-change of >1.20 or <0.83 and P-values < 0.05. Expression data were grouped together by hierarchical clustering
      according to the protein level. To annotate the sequences, information was extracted from UniProtKB/Swiss-Prot , 
      Kyoto Encyclopedia of Genes and Genomes (KEGG), and Gene Ontology (GO). GO and KEGG enrichment analyses were 
      carried out with the Fisher’s exact test, and FDR correction for multiple testing was also performed. GO terms 
      were grouped into three categories: biological process (BP), molecular function (MF), and cellular component (CC).
       Enriched GO and Kegg pathways were nominally statistically significant at the p&lt;0.05 level. Construction of 
       protein–protein interaction (PPI) networks were also conducted by using the STRING database with the cytoscape 
       software.</p> 
	<h1>
		<div>{1}. 拜谱服务产品<a name="lit2"></a></div>
	</h1>
	<div class="text-center">
      <figure class="figure">
           <a href="img/拜谱服务产品.png" class="thumbnail" target=_blank>
            <img src="img/拜谱服务产品.png" class="figure-img img-fluid d-block mx-auto" width=800 alt="无该结果"  />
            </a>
      </figure>
    </div>

	<h1>
		<div>{2}. 参考文献<a name="lit3"></a></div>
	</h1>
	<h6>[1]Cox, J. and M. Mann (2008). "MaxQuant enables high peptide identification rates, individualized p.p.b.-range
	 mass accuracies and proteome-wide protein quantification." Nat Biotechnol 26(12): 1367-1372.</h6>
	<h6>[2]Cox, J., N. Neuhauser, et al. (2011). "Andromeda: a peptide search engine integrated into the MaxQuant 
	environment." J Proteome Res 10(4): 1794-1805.</h6>
	<h6>[3]Luber CA, Cox J, Lauterbach H, Fancke B, Selbach M, Tschopp J, Akira S, Wiegand M, Hochrein H, 
	O'Keeffe M, Mann M. (2010) Quantitative proteomics reveals subset-specific viral recognition in dendritic
	 cells. Immunity 32, 279-89.</h6>
	<h6>[4]Waanders, LF; Chwalek K; Monetti M; Kumar C; Lammert E; Mann M; (2009) Quantitative proteomic analysis
	 of single pancreatic islets. PNAS 106, 18902-18907.</h6>
	<h6>[5]Ashburner M, Ball CA, et al. Gene ontology: tool for the unification of biology. The Gene Ontology 
	Consortium. Nat Genet. 2000; 25(1): 25-9.</h6>
	<h6>[6]Kanehisa M, Goto S, et al. KEGG for integration and interpretation of large-scale molecular data sets. 
	Nucleic Acids Res. 2012; 40(Database issue): D109-14.</h6>
	
	<h1>
		<div>{3}. 附录<a name="lit4"></a></div>
	</h1>
		<h4>|--1.DatabaseMatching/  搜库结果目录</h4>
		 <h5>|&nbsp;&nbsp;|--Protein.xlsx  蛋白鉴定表（表头说明见下表1）</h5>
		 <h5>|&nbsp;&nbsp;|--Peptide.xlsx  肽段鉴定表（表头说明见下表2）</h5>
		<h4>|-- 2.QualityControl/  数据质控目录</h4>
		 <h5>|&nbsp;&nbsp;|--peptide_error.pdf  肽段匹配误差分布图</h5>
		 <h5>|&nbsp;&nbsp;|--peptide_score.pdf  肽段得分分布柱状图</h5>
		 <h5>|&nbsp;&nbsp;|--peptide_length.pdf  肽段长度分布柱状图</h5>
		 <h5>|&nbsp;&nbsp;|--protein_MW&CI.pdf  蛋白分子量和等电点分布散点图</h5>
		 <h5>|&nbsp;&nbsp;|--protein_sequence_coverage.pdf  蛋白覆盖度分布环图</h5>
		 <h5>|&nbsp;&nbsp;|--protein_unique_peptide.pdf  unique肽段数量分布柱状图</h5>
		<h4>|-- 3.SampleAnalysis/  样本比较分析目录</h4>
		 <h5>|&nbsp;&nbsp;|--boxplot.pdf   样本表达分析箱线图</h5>
		 <h5>|&nbsp;&nbsp;|--violin.pdf   样本表达分析小提琴图</h5>
		 <h5>|&nbsp;&nbsp;|--Histogram.pdf   样本表达分析直方图</h5>
		 <h5>|&nbsp;&nbsp;|--Correlation_bubble.pdf   样本相关性分析气泡图</h5>
		 <h5>|&nbsp;&nbsp;|--Correlation_bcatter.pdf   样本相关性分析散点图</h5>
		 <h5>|&nbsp;&nbsp;|--PCA1.2.pdf, PCA1.2_style2.pdf, PCA1.2_style3 
		 不同风格的样本第一主成分和第二主成分的PCA分析得分图</h5>
         <h5>|&nbsp;&nbsp;|--PCA123.pdf 样本第一,第二，第三主成分的PCA分析得分图</h5>
		<h4>|-- 4.AllProtein/  全蛋白表达与功能注释结果目录</h4>
		  <h6>|&nbsp;&nbsp;|--summary.xlsx  全单白注释结果汇总表（表头说明见下表3）</h6>
		 <h5>|&nbsp;&nbsp;|--4.1GO/  GO term注释与富集分析结果目录（模式生物）</h5>
		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--GOlevel2CountBar.pdf   GO term注释统计柱状图（level2)</h6>
		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--pCountPoint.pdf  GO term富集气泡图（3大分支，横坐标p-value）</h6>
		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--pRFPoint.pdf  GO term富集气泡图（3大分支，横坐标RF）</h6>
		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--MF/BP/CC.EnrichedBar.pdf  GO term富集柱状图（单独分支</h6>
		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--MF/BP/CC.RichFactor.pdf  GO term top10富集气泡图（单独分支）</h6>
		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--MF/BP/CC.EnrichedSymbol.pdf  GO term top20富集气泡图（单独分支）</h6>
		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--MF/BP/CC.DAG.pdf  GO term富集有向无循环图（单独分支）</h6>
		 <h5>|&nbsp;&nbsp;|--4.2KEGG/  KEGG通路注释与富集分析结果目录（模式生物）</h5>
		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--Top10EnrichedBar.pdf     KEGG通路富集柱状图（top10）</h6>
		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--Top10EnrichSymbol.pdf     Top10 KEGG通路富集气泡图（top10）</h6>
		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--Top20EnrichedSymbol2.pdf       Top20 KEGG通路富集气泡图（top20）</h6>
		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--KEGG.Sig.Bar.pdf  Top30 KEGG通路富集柱状图（颜色归属level1）</h6>
		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--KEGG.Sig.Bar2.pdf      Top30 KEGG通路富集柱状图（颜色归属level2）</h6>
		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--kegg.sig.pdf  KEGG通路富集top20气泡图（显示level 1分类）</h6>
		 <h5>|&nbsp;&nbsp;|--4.3subcellular_localization/  亚细胞定位（CC）分析结果目录</h5>
		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--CellLocationMerge.png    亚细胞定位统计结果（柱状图、点图、饼图、环图）</h6>
		<h4>|-- 5.DEPStatistics/  表达差异分析结果目录</h4>
		 <h5>|&nbsp;&nbsp;|--DEP. xlsx    差异表达蛋白统计表（表头说明见下表5）</h5>
         <h5>|&nbsp;&nbsp;|--Venn/   Venn分析结果（至少2组比较)</h5>
		 <h5>|&nbsp;&nbsp;|--Venn.pdf          各组蛋白鉴定数量韦恩图</h5>
		 <h5>|&nbsp;&nbsp;|--Venn.xlsx          韦恩分析统计表（表头说明见下表6）</h5>
         <h5>|&nbsp;&nbsp;|--*.vs.*  各比对组结果</h5>
		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.volcano.pdf  各比对组差异蛋白可视化火山图</h6>
		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.FoldChange.scatter.pdf  各比对组蛋白FC可视化散点图</h6>
		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.pvalue.scatter.pdf  各比对组蛋白p-value可视化散点图</h6>
		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.updown.bar.pdf  各比对组上下调蛋白统计柱状图</h6>
		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.DEP.Donut.pdf  各比对组差异蛋白统计环图</h6>
		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.heatmap.pdf  各比对组差异蛋白聚类热图（两组）</h6>
		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.heatmap-Multi.pdf  各比对组差异蛋白聚类热图（至少三组）</h6>
		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.heatmap.xlsx  各比对组差异蛋白聚类统计表（表头说明见下表7）</h6>
		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.K-value trend.pdf  各比对组差异蛋白k值聚类趋势图（至少三组</h6>
		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.K-value.xlsx  各比对组差异蛋白k值聚类统计表（表头说明见下表8）</h6>
		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.BarComp/Boxplot/Jitter/ScatterSE/ 
		 各比对组差异蛋白表达分析结果目录：比较柱状图/箱线图/抖图/散点图（两组)</h6>
		<h4>|-- 6.DEPFunction/  差异蛋白功能分析结果目录</h4>
         <h5>|&nbsp;&nbsp;|--*.vs.*      比较组的结果</h5>
		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.[all/up/down]_summary.xlsx 
		 all/up/down差异蛋白注释汇总表（表头说明见下表9）</h6>
		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.[all/up/down]_GO_Ident/ 
		 以全部鉴定到的蛋白为背景参照的差异蛋白GO分析结果目录</h6>
         <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.[all/up/down]_GO_Species/ 
         以该物种全部基因为背景参照的差异蛋白GO分析结果目录</h6>
         <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.[all/up/down]_GO_Level2/  第二层次GO功能的统计结果</h6>
         <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.[all/up/down]_subcellular_localization/ 
         经典的亚细胞器定位的统计结果</h6>
         <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.[all/up/down]_KEGG  差异蛋白KEGG分析结果目录</h6>
		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.all/up/down_PFIN/  all/up/down差异蛋白功能互作（PFI）分析目录</h6>
		<h4>|-- 7.Anova/  多组比较分析结果（如果存在多组比较分析）</h4>
		  <h5>|&nbsp;&nbsp;|--[比较组的名称,由组别加'_'组成<h5>
		   <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--Multi_group_expression         多组表达单个蛋白绘图</h6>
		   <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--kmeans         k-means分析结果</h6>
		   <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--GOEnrich_Identified         以鉴定到的全蛋白为背景参照的GO富集分析结果</h6>
		   <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--GOEnrich_Species         以该物种全部基因为背景参照的GO富集分析结果</h6>
		   <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--GOLevel2         第二层GO功能的统计结果</h6>
		   <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--GO_subcellular_localization         蛋白质在经典的亚细胞器定位的统计结果</h6>
		   <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--KEGG         KEGG通路分析结果</h6>
	
    <div class="gohome">回到顶部</div>
  <script src="js/back.js"></script>
""".format(self.index, self.index + 1, self.index + 2, self.index + 3))

        self.print_bottom(ero)
        url='literature.html'
        self.lefthandle.write("""<div class="level1"><a href="{}#{}" target="rightframe">{}</div>""".format(url, 'lit1',
                                                                                                            str(
                                                                                                                self.index) + '. ' + "英文写作参考"))
        self.lefthandle.write("""<div class="level1"><a href="{}#{}" target="rightframe">{}</div>""".format(url, 'lit2',
                                                                                                            str(
                                                                                                                self.index+1) + '. ' + "拜谱服务产品"))
        self.lefthandle.write("""<div class="level1"><a href="{}#{}" target="rightframe">{}</div>""".format(url, 'lit3',
                                                                                                            str(
                                                                                                                self.index + 2) + '. ' + "参考文献"))
        self.lefthandle.write("""<div class="level1"><a href="{}#{}" target="rightframe">{}</div>""".format(url, 'lit4',
                                                                                                            str(
                                                                                                                self.index + 3) + '. ' + "附录"))
    def print_lfq_english(self):
        ero = open(self.pn + '/html/literature.html', 'w')
        self.print_head(ero)
        ero.write("""
        	<div id="english">
        		<h1>
        			<div>{0}. 英文写作参考<a name="lit1"></a></div>
        		</h1>
        	</div>	
        	<div class="alert alert-warning">注：下文仅作为英文模板参考，具体细节要以中文报告为准！</div>
        	<h2>Materials and Methods</h2>
        	<h2>Reagents</h2>
        	<p>Ammonium bicarbonate, dithiothreitol (DTT), iodoacetamide (IAA), and sodium carbonate were purchased 
        	from Sigma-Aldrich (St. Louis, MO). Urea and Sodium dodecyl sulfate (SDS) were purchased from Bio-Rad 
        	(Hercules, CA). Acetonitrile and water for nano-LC−MS/MS were purchased from J. T. Baker (Phillipsburg, NJ).
        	 Trypsin was purchased from Promega (Madison, WI). All other chemical reagent were purchased with analytical
        	  grade.</p>
        	<h2>Sample preparation</h2>
        	<p>Protein was extracted from tissue samples using SDT lysis buffer (4% SDS, 100 mM DTT, 100 mM Tris-HCl 
        	pH 8.0). Sample were boiled for 5 min and further ultrasonicated and boiling again for another 5 min. 
        	Undissolved cellular debris were removed by centrifugation at 16 000g for 15 min. The supernatant was 
        	collected and quantified with a BCA Protein Assay Kit (Bio-Rad, USA).</p>

        	<h2>Protein Digestion</h2>
        	<p>Protein (200 μg for each sample) digestion was performed with FASP method described by Wisniewski, 
        	Zougman et al[6]. Briefly, the detergent, DTT and IAA in UA buffer was added to block reduced cysteine. 
        	Finally, the protein suspension was digested with trypsin (Promega) at ratio 50:1 overnight at 37 °C. 
        	The peptide were collected by centrifugation at 16 000g for 15 min. The peptide was desalted with C18 
        	StageTip for further LC-MS analysis. The concentrations of peptides were determined with OD280 by Nanodrop 
        	One device. </p>
        	<h2>LC-MS/MS Analysis</h2>
        	<p>LC-MS/MS were performed on a Q Exactive Plus mass spectrometer coupled with Easy 1200 nLC (Thermo Fisher
        	 Scientific). Peptide was first loaded to a trap column (100μm*20mm, 5μm，C18，Dr. Maisch GmbH, Ammerbuch, 
        	 Germany) in buffer A (0.1% Formic acid in water). Reverse-phase high-performance liquid chromatography 
        	 (RP-HPLC) separation was performed with the EASY-nLC system (Thermo Fisher Scientific, Bremen, Germany)
        	  using a self-packed column (75 μm × 150 mm; 3 μm ReproSil-Pur C18 beads, 120 Å, Dr. Maisch GmbH, 
        	  Ammerbuch, Germany) at a flow rate of 300 nL/min. The RP−HPLC mobile phase A was 0.1% formic acid in 
        	  water, and B was 0.1% formic acid in 95% acetonitrile. Peptide were eluted over 120 min with a linear 
        	  gradient of buffer B. MS data was acquired using a data-dependent top20 method dynamically choosing the 
        	  most abundant precursor ions from the survey scan (300–1800 m/z) for HCD fragmentation. The instrument 
        	  was run with peptide recognition mode enabled. A lock mass of 445.120025 Da was used as internal standard 
        	  for mass calibration. The full MS scans were acquired at a resolution of 70,000 at m/z 200, and 17,500 at 
        	  m/z 200 for MS/MS scan. The maximum injection time was set to for 50 ms for MS and 50 ms for MS/MS. 
        	  Normalized collision energy was 27 and the isolation window was set to 1.6 Th. Dynamic exclusion 
        	  duration was 60 s.</p>	
        	<h2>Sequence Database Searching and Data Analysis</h2>
            <p>The MS data were analyzed using MaxQuant software version 1.6.0.16. MS data were searched against the
             UniProtKB Rattus norvegicus database (36080 total entries, downloaded 08/14/2018). (根据数据库情况修改)。
             The trypsin was seleted as digestion enzyme. The maximal two missed cleavage sites and the mass tolerance 
             of 4.5 ppm for precursor ions and 20 ppm for fragment ions were defined for database search. 
             Carbamidomethylation of cysteines was defined as fixed modification, while acetylation of protein 
             N-terminal, oxidation of Methionine was set as variable modifications for database searching. The database
              search results were filtered and exported with <1% false discovery rate (FDR) at peptide-spectrum-matched
               level, and protein level, respectively Label-free quantification was carried out in MaxQuant using 
               intensity determination and normalization algorithm as previously described [7, 8, 9]. 
               The “LFQ intensity” of each protein in different samples was calculated as the best estimate, satisfying
                all of the pairwise peptide comparisons, and this LFQ intensity was almost on the same scale of the 
                summed-up peptide intensities. The quantitative protein ratios were weighted and normalized by the 
                median ratio in Maxquant software. Only proteins with fold change ≥1.5-fold and a p-value <0.05 were 
                considered for significantly differential expressions.</p>	
            <h2>Bioinformatics analysis</h2>
            <p>Analyses of bioinformatics data were carried out with Perseus software [10], Microsoft Excel and R 
            statistical computing software. Hierarchical clustering analysis was performed with the pheatmap package, 
            which is based on the open-source statistical language R25, using Euclidean distance as the distance metric 
            and complete method as the agglomeration method. To annotate the sequences, information was extracted from 
            UniProtKB/Swiss-Prot [11], Kyoto Encyclopedia of Genes and Genomes (KEGG) [12], and Gene Ontology (GO) 
            [13]. GO and KEGG enrichment analyses were carried out with the Fisher’s exact test, and FDR correction for 
            multiple testing was also performed. GO terms were grouped into three categories: biological process (BP), 
            molecular function (MF), and cellular component (CC) [13]. Enriched GO and Kegg pathways were nominally 
            statistically significant at the p&lt;0.05 level. Construction of protein–protein interaction (PPI) networks 
            were also conducted by using the STRING database with the Cytoscape software [14].</p> 
        	<h1>{1}. 拜谱服务产品<a name="lit2"></a></h1>
        	<div class="text-center">
              <figure class="figure">
                   <a href="img/拜谱服务产品.png" class="thumbnail" target=_blank>
                    <img src="img/拜谱服务产品.png" class="figure-img img-fluid d-block mx-auto" width=800 alt="无该结果"  />
                    </a>
              </figure>
            </div>

        	<h1>{2}. 参考文献<a name="lit3"></a></h1>
        	<h6>[1]Cox, J. and M. Mann (2008). "MaxQuant enables high peptide identification rates, individualized p.p.b.-range
        	 mass accuracies and proteome-wide protein quantification." Nat Biotechnol 26(12): 1367-1372.</h6>
        	<h6>[2]Cox, J., N. Neuhauser, et al. (2011). "Andromeda: a peptide search engine integrated into the MaxQuant 
        	environment." J Proteome Res 10(4): 1794-1805.</h6>
        	<h6>[3]	Cox, J., M. Y. Hein, et al. (2014). "Accurate proteome-wide label-free quantification by delayed 
        	normalization and maximal peptide ratio extraction, termed MaxLFQ." Mol Cell Proteomics 13(9): 2513-2526.</h6>
        	<h6>[4]Ashburner M, Ball CA, et al. Gene ontology: tool for the unification of biology. The Gene Ontology 
        	Consortium. Nat Genet. 2000; 25(1): 25-9.</h6>
        	<h6>[5]Kanehisa M, Goto S, et al. KEGG for integration and interpretation of large-scale molecular data sets. 
        	Nucleic Acids Res. 2012; 40(Database issue): D109-14.</h6>
        	<h6>[6]Wisniewski, J.R., et al., Universal sample preparation method for proteome analysis. Nat Methods, 
        	2009. 6(5): p. 359-62.</h6>
        	<h6>[7]	Schwanhausser, B., et al., Global quantification of mammalian gene expression control. Nature, 
        	2011. 473(7347): p. 337-42.</h6>
<h6>[8]	Luber, C.A., et al., Quantitative proteomics reveals subset-specific viral recognition in dendritic cells. 
Immunity, 2010. 32(2): p. 279-89.</h6>
<h6>[9]	Cox, J., M. Y. Hein, et al. (2014). "Accurate proteome-wide label-free quantification by delayed normalization 
and maximal peptide ratio extraction, termed MaxLFQ." Mol Cell Proteomics 13(9): 2513-2526.</h6>
<h6>[10] Tyanova S, Temu T, Sinitcyn P. The Perseus computational platform for comprehensive analysis of 
(prote)omics data. 2016;13(9):731-740.</h6>
<h6>[11] Boutet E, Lieberherr D, Tognolli M, et al. UniProtKB/Swiss-Prot, the Manually Annotated Section of the 
UniProt KnowledgeBase: How to Use the Entry View. Methods in molecular biology (Clifton, NJ) 2016;1374:23-54.</h6>
<h6>[12] Kanehisa M, Goto S, Sato Y, Furumichi M, Tanabe M. KEGG for integration and interpretation of large-scale 
molecular data sets. Nucleic acids research 2012;40(Database issue):D109-114.</h6>
<h6>[13] Ashburner M, Ball CA, Blake JA, et al. Gene ontology: tool for the unification of biology. The Gene Ontology 
Consortium. Nature genetics 2000;25(1):25-29.</h6>
<h6>[14] Kohl M, Wiese S, Warscheid B. Cytoscape: software for visualization and analysis of biological networks. 
Methods in molecular biology (Clifton, NJ) 2011;696:291-303.</h6>
        	<h1>{3}. 附录<a name="lit4"></a></h1>
        		<h4>|--1.DatabaseMatching/  搜库结果目录</h4>
        		 <h5>|&nbsp;&nbsp;|--Protein.xlsx  蛋白鉴定表（表头说明见下表1）</h5>
        		 <h5>|&nbsp;&nbsp;|--Peptide.xlsx  肽段鉴定表（表头说明见下表2）</h5>
        		<h4>|-- 2.QualityControl/  数据质控目录</h4>
        		 <h5>|&nbsp;&nbsp;|--peptide_error.pdf  肽段匹配误差分布图</h5>
        		 <h5>|&nbsp;&nbsp;|--peptide_score.pdf  肽段得分分布柱状图</h5>
        		 <h5>|&nbsp;&nbsp;|--peptide_length.pdf  肽段长度分布柱状图</h5>
        		 <h5>|&nbsp;&nbsp;|--protein_MW&CI.pdf  蛋白分子量和等电点分布散点图</h5>
        		 <h5>|&nbsp;&nbsp;|--protein_sequence_coverage.pdf  蛋白覆盖度分布环图</h5>
        		 <h5>|&nbsp;&nbsp;|--protein_unique_peptide.pdf  unique肽段数量分布柱状图</h5>
        		<h4>|-- 3.SampleAnalysis/  样本比较分析目录</h4>
        		 <h5>|&nbsp;&nbsp;|--boxplot.pdf   样本表达分析箱线图</h5>
        		 <h5>|&nbsp;&nbsp;|--violin.pdf   样本表达分析小提琴图</h5>
        		 <h5>|&nbsp;&nbsp;|--Histogram.pdf   样本表达分析直方图</h5>
        		 <h5>|&nbsp;&nbsp;|--Correlation_bubble.pdf   样本相关性分析气泡图</h5>
        		 <h5>|&nbsp;&nbsp;|--Correlation_bcatter.pdf   样本相关性分析散点图</h5>
        		 <h5>|&nbsp;&nbsp;|--PCA1.2.pdf, PCA1.2_style2.pdf, PCA1.2_style3   不同风格的样本第一主成分和第二主成分的PCA分析得分图</h5>
                 <h5>|&nbsp;&nbsp;|--PCA123.pdf 样本第一,第二，第三主成分的PCA分析得分图</h5>
        		<h4>|-- 4.AllProtein/  全蛋白表达与功能注释结果目录</h4>
        		  <h6>|&nbsp;&nbsp;|--summary.xlsx  全单白注释结果汇总表（表头说明见下表3）</h6>
        		 <h5>|&nbsp;&nbsp;|--4.1GO/  GO term注释与富集分析结果目录（模式生物）</h5>
        		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--GOlevel2CountBar.pdf   GO term注释统计柱状图（level2)</h6>
        		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--pCountPoint.pdf  GO term富集气泡图（3大分支，横坐标p-value）</h6>
        		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--pRFPoint.pdf  GO term富集气泡图（3大分支，横坐标RF）</h6>
        		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--MF/BP/CC.EnrichedBar.pdf  GO term富集柱状图（单独分支</h6>
        		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--MF/BP/CC.RichFactor.pdf  GO term top10富集气泡图（单独分支）</h6>
        		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--MF/BP/CC.EnrichedSymbol.pdf  GO term top20富集气泡图（单独分支）</h6>
        		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--MF/BP/CC.DAG.pdf  GO term富集有向无循环图（单独分支）</h6>
        		 <h5>|&nbsp;&nbsp;|--4.2KEGG/  KEGG通路注释与富集分析结果目录（模式生物）</h5>
        		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--Top10EnrichedBar.pdf     KEGG通路富集柱状图（top10）</h6>
        		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--Top10EnrichSymbol.pdf     Top10 KEGG通路富集气泡图（top10）</h6>
        		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--Top20EnrichedSymbol2.pdf       Top20 KEGG通路富集气泡图（top20）</h6>
        		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--KEGG.Sig.Bar.pdf  Top30 KEGG通路富集柱状图（颜色归属level1）</h6>
        		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--KEGG.Sig.Bar2.pdf      Top30 KEGG通路富集柱状图（颜色归属level2）</h6>
        		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--kegg.sig.pdf  KEGG通路富集top20气泡图（显示level 1分类）</h6>
        		 <h5>|&nbsp;&nbsp;|--4.3subcellular_localization/  亚细胞定位（CC）分析结果目录</h5>
        		  <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--CellLocationMerge.png    亚细胞定位统计结果（柱状图、点图、饼图、环图）</h6>
        		<h4>|-- 5.DEPStatistics/  表达差异分析结果目录</h4>
        		 <h5>|&nbsp;&nbsp;|--DEP. xlsx    差异表达蛋白统计表（表头说明见下表5）</h5>
                 <h5>|&nbsp;&nbsp;|--Venn/   Venn分析结果（至少2组比较)</h5>
        		 <h5>|&nbsp;&nbsp;|--Venn.pdf          各组蛋白鉴定数量韦恩图</h5>
        		 <h5>|&nbsp;&nbsp;|--Venn.xlsx          韦恩分析统计表（表头说明见下表6）</h5>
                 <h5>|&nbsp;&nbsp;|--*.vs.*  各比对组结果</h5>
        		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.volcano.pdf  各比对组差异蛋白可视化火山图</h6>
        		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.FoldChange.scatter.pdf  各比对组蛋白FC可视化散点图</h6>
        		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.pvalue.scatter.pdf  各比对组蛋白p-value可视化散点图</h6>
        		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.updown.bar.pdf  各比对组上下调蛋白统计柱状图</h6>
        		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.DEP.Donut.pdf  各比对组差异蛋白统计环图</h6>
        		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.heatmap.pdf  各比对组差异蛋白聚类热图（两组）</h6>
        		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.heatmap-Multi.pdf  各比对组差异蛋白聚类热图（至少三组）</h6>
        		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.heatmap.xlsx  各比对组差异蛋白聚类统计表（表头说明见下表7）</h6>
        		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.K-value trend.pdf  各比对组差异蛋白k值聚类趋势图（至少三组</h6>
        		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.K-value.xlsx  各比对组差异蛋白k值聚类统计表（表头说明见下表8）</h6>
        		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.BarComp/BarSE/Boxplot/Jitter/ScatterSE/    各比对组差异蛋白表达分析结果目录：柱状图/箱线图/抖图/散点图（两组</h6>
        		<h4>|-- 6.DEPFunction/  差异蛋白功能分析结果目录</h4>
                 <h5>|&nbsp;&nbsp;|--*.vs.*      比较组的结果</h5>
        		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.[all/up/down]_summary.xlsx      all/up/down差异蛋白注释汇总表（表头说明见下表9）</h6>
        		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.[all/up/down]_GO_Ident/  以全部鉴定到的蛋白为背景参照的差异蛋白GO分析结果目录</h6>
                 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.[all/up/down]_GO_Species/  以该物种全部基因为背景参照的差异蛋白GO分析结果目录</h6>
                 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.[all/up/down]_GO_Level2/  第二层次GO功能的统计结果</h6>
                 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.[all/up/down]_subcellular_localization/  经典的亚细胞器定位的统计结果</h6>
                 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.[all/up/down]_KEGG  差异蛋白KEGG分析结果目录</h6>

        		 <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--*.vs.*.all/up/down_PFIN/  all/up/down差异蛋白功能互作（PFI）分析目录</h6>
        	<h4>|-- 7.Anova/  多组比较分析结果（如果存在多组比较分析）</h4>
		  <h5>|&nbsp;&nbsp;|--[比较组的名称,由组别加'_'组成<h5>
		   <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--Multi_group_expression         多组表达单个蛋白绘图</h6>
		   <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--kmeans         k-means分析结果</h6>
		   <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--GOEnrich_Identified         以鉴定到的全蛋白为背景参照的GO富集分析结果</h6>
		   <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--GOEnrich_Species         以该物种全部基因为背景参照的GO富集分析结果</h6>
		   <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--GOLevel2         第二层GO功能的统计结果</h6>
		   <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--GO_subcellular_localization         蛋白质在经典的亚细胞器定位的统计结果</h6>
		   <h6>|&nbsp;&nbsp;|&nbsp;&nbsp;|--KEGG         KEGG通路分析结果</h6>
            <div class="gohome">回到顶部</div>
          <script src="js/back.js"></script>
        """.format(self.index, self.index + 1, self.index + 2, self.index + 3))

        self.print_bottom(ero)
        url = 'literature.html'
        self.lefthandle.write("""<div class="level1"><a href="{}#{}" target="rightframe">{}</div>""".format(url, 'lit1',
                                                                                                            str(
                                                                                                                self.index) + '. ' + "英文写作参考"))
        self.lefthandle.write("""<div class="level1"><a href="{}#{}" target="rightframe">{}</div>""".format(url, 'lit2',
                                                                                                            str(
                                                                                                                self.index + 1) + '. ' + "拜谱服务产品"))
        self.lefthandle.write("""<div class="level1"><a href="{}#{}" target="rightframe">{}</div>""".format(url, 'lit3',
                                                                                                            str(
                                                                                                                self.index + 2) + '. ' + "参考文献"))
        self.lefthandle.write("""<div class="level1"><a href="{}#{}" target="rightframe">{}</div>""".format(url, 'lit4',
                                                                                                            str(
                                                                                                                self.index + 3) + '. ' + "附录"))

    def close(self):
        self.lefthandle.write("</div></body></html>")
        self.lefthandle.close()

    def print_firstpage(self):
        tp = open(self.pn+'/html/firstpage.html', 'w')
        self.print_head_left(tp)
        tp.write("""<div><a name="fengmian"></a><img src="img/fengmian.png" width = 100% /><a name="fengmian"></a></div>
        
        """)
        self.print_h1(tp, "主要分析内容", "firstpage.html", "firstpage")
        tp.write("""<div class="text-center"><img src="img/anal.png" width=800 /></div>""")
        self.print_bottom(tp)
        tp.close()

        self.index = self.index + 1

if __name__ == '__main__':
    bin = get_absdir()
    # lisense check
    licensefile = bin + '/License.dat'
    hylicense.eda.license_check(licensefile, 'prot')

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='ifold', default='tmpinfold', help='挑选好的结果的目录')
    parser.add_argument('-ia', action='store', dest='pn', default='tmp',
                        help='总的结果目录')
    parser.add_argument('-t', action='store', dest='type', default='TMT',
                        help='报告类型,TMT|LFQ|DIA,default is TMT')
    p = parser.parse_args()

    if not os.path.exists(p.ifold) or not os.path.exists(p.pn):
        parser.print_help()
        sys.exit()



    report = protreport(p.ifold, p.ifold + "/report.html", p.pn)

    #拷贝表格说明
    if os.path.exists(bin + "/module/table_readme.pdf"):
        cpfile(bin + "/module/table_readme.pdf", report.pn)

    #拷贝制裁支持文件
    if p.type == 'TMT':
        if not os.path.exists(report.pn+'/html/CSS'):
            shutil.copytree(bin+"/module/tmtreport/CSS", report.pn+'/html/CSS')
        if not os.path.exists(report.pn+'/html/img'):
            shutil.copytree(bin+"/module/tmtreport/img", report.pn+'/html/img')
        if not os.path.exists(report.pn+'/html/js'):
            shutil.copytree(bin+"/module/tmtreport/js", report.pn+'/html/js')
    elif p.type == 'LFQ':
        if not os.path.exists(report.pn+'/html/CSS'):
            shutil.copytree(bin+"/module/lfqreport/CSS", report.pn+'/html/CSS')
        if not os.path.exists(report.pn+'/html/img'):
            shutil.copytree(bin+"/module/lfqreport/img", report.pn+'/html/img')
        if not os.path.exists(report.pn+'/html/js'):
            shutil.copytree(bin+"/module/lfqreport/js", report.pn+'/html/js')
    elif p.type == 'DIA':
        if not os.path.exists(report.pn+'/html/CSS'):
            shutil.copytree(bin+"/module/lfqreport/CSS", report.pn+'/html/CSS')
        if not os.path.exists(report.pn+'/html/img'):
            shutil.copytree(bin+"/module/lfqreport/img", report.pn+'/html/img')
        if not os.path.exists(report.pn+'/html/js'):
            shutil.copytree(bin+"/module/lfqreport/js", report.pn+'/html/js')
        os.remove(report.pn+'/html/img/fengmian.png')
        shutil.copy(bin+"/module/fengmian.png", report.pn+'/html/img')
    else:
        pass

    report.print_firstpage()
    report.print_all_protein_function_analysis()
    report.print_dep_protein_function_analysis()
    report.print_multi_group_comparison_analysis(report.index)

    report.print_advance_html()
    report.print_personal_analysis()
    report.print_qc_analysis()
    report.print_sample_compare_analysis()
    if p.type == 'TMT':
        report.print_english_reference()
    else:
        report.print_lfq_english()

    report.close()