# -*- coding: utf-8 -*-
import os
import re
import sys
import time
from copy import copy
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from scGBM_ui import Ui_scGBM
from PyQt5 import Qt
from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import QCompleter, QMessageBox, QFileDialog, QCheckBox, QHBoxLayout, QWidget, QTableWidgetItem, QAbstractItemView
# from platform_childthread import Message
import pandas as pd
import numpy as np
import seaborn as sns
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import zipfile
import os
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import pearsonr
from scipy.stats import mannwhitneyu

matplotlib.use('Qt5Agg')
gui_continue = QtWidgets.QApplication.processEvents

class myGUI(Ui_scGBM):
    def __init__(self):
        "Initialize"
        self.threadpool = {}
        self.filepath = '/data/activate_data/majiayan/scGBM_data/gbm.pkl'
        self.tcgapath = '/data/activate_data/majiayan/scGBM_data/TCGA_Survival.txt'
        self.cggapath = '/data/activate_data/majiayan/scGBM_data/CGGA_Survival.txt'
        self.data = {}
        self.tcgadata = []
        self.cggadata = []
        self.genelist = []
        self.plotdict = {}
        self.statdict = {}
        self.running_report = ''
        self.message = ''
        self.lines = []
    # UI
    def set_default(self):
        "Set default value"
        # 设置默认值
        self.pointsize_box.setRange(1, 10);
        self.pointsize_slider.setRange(1, 10);
        self.pointsize_slider.setValue(1)
        self.pointsize_box.setValue(1)
        self.read_data_start()
        self.SubtypeTable.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.FiguresTab.setCurrentIndex(0)

    def setup_signal(self):
        """Connect signal/slot"""
        # 主窗口上半部分
        #self.filter_1.currentTextChanged.connect(self.filter_2_refresh)  # 刷新类型选择
        self.Run.clicked.connect(self.preview)  # 预览
        self.Save.clicked.connect(self.save_fig)  # 保存图片
        self.pointsize_slider.valueChanged.connect(lambda:self._splider_change())#滑块的connect
        self.pointsize_box.valueChanged.connect(lambda:self._spinbox_change())#微调框的connect
        self.SampleNameBox.currentIndexChanged.connect(self.completer_init)
        self.SampleNameBox.currentIndexChanged.connect(self.setTable)
        self.SampleNameBox.currentIndexChanged.connect(self.setGene)

    def _splider_change(self):
        self.pointsize_box.setValue(self.pointsize_slider.value())
      
    def _spinbox_change(self):
        self.pointsize_slider.setValue(self.pointsize_box.value())

    def setTable(self):
        self.lines=[]
        row = 0
        self.SubtypeTable.setRowCount(row+1)
        ck = QCheckBox()
        ck.setChecked(True)
        h = QHBoxLayout()
        h.setAlignment(QtCore.Qt.AlignCenter)
        h.addWidget(ck)
        w = QWidget()
        w.setLayout(h)
        totalcount = 0
        for i in list(self.data[self.SampleNameBox.currentText()].keys()):
            totalcount = totalcount + len(self.data[self.SampleNameBox.currentText()][i]['cellType'].tolist())
        cellcount = str(totalcount)
        #print(cellcount)
        self.SubtypeTable.setCellWidget(row,0,w)
        self.SubtypeTable.setItem(row,1,QTableWidgetItem('all'))
        self.SubtypeTable.setItem(row,2,QTableWidgetItem(cellcount))
        self.lines.append([ck,'all',cellcount])
        row += 1
        for celltype in list(self.data[self.SampleNameBox.currentText()].keys()):
            #print(celltype)
            self.SubtypeTable.setRowCount(row+1)
            ck = QCheckBox()
            h = QHBoxLayout()
            h.setAlignment(QtCore.Qt.AlignCenter)
            h.addWidget(ck)
            w = QWidget()
            w.setLayout(h)
            cellcount = str(len(self.data[self.SampleNameBox.currentText()][celltype]['cellType'].tolist()))
            #print(cellcount)
            self.SubtypeTable.setCellWidget(row,0,w)
            self.SubtypeTable.setItem(row,1,QTableWidgetItem(celltype))
            self.SubtypeTable.setItem(row,2,QTableWidgetItem(cellcount))
            row += 1
            self.lines.append([ck,celltype,cellcount])
        header = self.SubtypeTable.horizontalHeader()       
        header.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
        header.setSectionResizeMode(2, QtWidgets.QHeaderView.ResizeToContents)

    def backfeed_update(self, text):
        """backfeed widget update"""
        self.running_report += text
        self.running_report += '\n'
        self.backfeed.setText(self.running_report)
        self.backfeed.verticalScrollBar().setValue(self.backfeed.verticalScrollBar().maximum())

    def completer_init(self):
        """Initialize a QCompleter"""
        self.completer = QCompleter(sorted(self.data[self.SampleNameBox.currentText()][list(self.data[self.SampleNameBox.currentText()].keys())[0]].columns[3:]))
        self.completer.setCompletionMode(QCompleter.UnfilteredPopupCompletion)
        self.GeneName_Line.setCompleter(self.completer)
        self.completer2 = QCompleter(sorted(self.data[self.SampleNameBox.currentText()][list(self.data[self.SampleNameBox.currentText()].keys())[0]].columns[3:]))
        self.completer2.setCompletionMode(QCompleter.UnfilteredPopupCompletion)
        self.GeneName_Line_2.setCompleter(self.completer2)

    def type_gene(self):
        """Type gene in lineedit, automatic completion with self.completer"""
        if self.input_gene.text() not in self.genelist and self.input_gene.text() in self.data.var_names:
            self.genelist.insert(0, self.input_gene.text())
            self.input_gene_multi.setText('\n'.join(self.genelist))
            self.genelist_refresh()

    def genelist_refresh(self):
        """Input genes in textedit"""
        self.genelist = self.input_gene_multi.toPlainText().split('\n')
        while '' in self.genelist:
            self.genelist.remove('')
        self.genelist = [gene for gene in self.genelist if gene in self.data.var_names]

    def import_genes(self):
        filepath = QFileDialog.getOpenFileName(self.centralwidget, '选择文件', os.getcwd(), 'Text(*.txt)')[0]
        genelist = [gene.strip() for gene in open(filepath)]
        self.import_genelist.setText('\n'.join(genelist))
        self.genelist_refresh()

    def setGene(self):
        self.GeneName_Line.setText(self.data[self.SampleNameBox.currentText()][list(self.data[self.SampleNameBox.currentText()].keys())[0]].columns[3])
        self.GeneName_Line_2.setText('')

    # backend
    def read_data_start(self):
        """Start a QThread for read data"""
        #self.backfeed_update('Data loading')
        from scGBM_childthread import Loading
        loading_thread = Loading(self.filepath,self.tcgapath,self.cggapath)
        loading_thread.signal.connect(self.read_data_end)
        self.threadpool['Loading'] = loading_thread
        loading_thread.start()

    def read_data_end(self):
        """1. Assigned Qthread.data to self.data
        2. Fill NA
        3. Initialize completer"""
        loading = self.threadpool['Loading']
        self.data = loading.data
        self.tcgadata = loading.tcgadata
        self.cggadata = loading.cggadata
        #print(self.data.keys())
        self.SampleNameBox.addItems(self.data.keys())
        self.completer_init()
        # delete finished thread
        del self.threadpool['Loading']

    def filter_2_refresh(self):
        """Refresh filter2 when filter1 been choosed"""
        if self.filter_1.currentText() == '无':
            self.comboBox_filter2.setCurrentText('无')
        else:
            filter2 = list(set(self.data.obs[self.filter_1.currentText()]))
            print(filter2)
            filter2.sort()
            filter2.insert(0, '无')
            self.comboBox_filter2.addItems(filter2)
            self.comboBox_filter2.setCurrentText('无')

    def preview(self):
        """Preview plot
        1. violin plot
        2. bubble plot
        3. survival curve
        4. cox regrssion
        5. single cell series(high expression gene,)"""
        plt.close()
        subtypes = []
        for line in self.lines:
            if line[0].isChecked():
                subtypes.append(line[1])
        pattern = re.compile(r'\d+')

        print('umap')
        fig = self.drawumap(subtypes)
        width, height = pattern.findall(fig.__str__())
        canvas = FigureCanvas(fig)
        canvas.setMinimumSize(int(width), int(height))
        self.showumap.setWidget(canvas)

        print('violin')
        fig = self.drawViolin(subtypes)
        width, height = pattern.findall(fig.__str__())
        canvas = FigureCanvas(fig)
        canvas.setMinimumSize(int(width), int(height))
        self.showviolin.setWidget(canvas)

        print('bar')
        fig = self.drawBar(subtypes)
        width, height = pattern.findall(fig.__str__())
        canvas = FigureCanvas(fig)
        canvas.setMinimumSize(int(width), int(height))
        self.showbar.setWidget(canvas) 

        print('subumap')
        fig = self.drawsubumap(subtypes)
        width, height = pattern.findall(fig.__str__())
        canvas = FigureCanvas(fig)
        canvas.setMinimumSize(int(width), int(height))
        self.showsubumap.setWidget(canvas)    

        print('tcgasurv')
        fig = self.drawTCGASurvival()
        width, height = pattern.findall(fig.__str__())
        canvas = FigureCanvas(fig)
        canvas.setMinimumSize(int(width), int(height))
        self.showtcgasurvival.setWidget(canvas)

        print('cggasurv')
        fig = self.drawCGGASurvival()
        width, height = pattern.findall(fig.__str__())
        canvas = FigureCanvas(fig)
        canvas.setMinimumSize(int(width), int(height))
        self.showcggasurvival.setWidget(canvas)

        print('corr')
        fig = self.drawCorr(subtypes)
        width, height = pattern.findall(fig.__str__())
        canvas = FigureCanvas(fig)
        canvas.setMinimumSize(int(width), int(height))
        self.showgenecorr.setWidget(canvas)   

        print('de')
        fig = self.calcDE(subtypes)
        width, height = pattern.findall(fig.__str__())
        canvas = FigureCanvas(fig)
        canvas.setMinimumSize(int(width), int(height))
        self.showde.setWidget(canvas)   

    def calcDE(self,subtypes):

        fig, axes = plt.subplots(1,1)

        if 'all' in subtypes:
            celllist = list(self.data[self.SampleNameBox.currentText()])
        else:
            celllist = subtypes
        gene = self.GeneName_Line.text()

        corr = np.array([[1.0]*len(celllist)]*len(celllist))
        for i in range(len(celllist)-1):
            for j in range(i,len(celllist)):
                if i == j:
                    continue
                if sum(self.data[self.SampleNameBox.currentText()][celllist[i]][gene].tolist()) == 0 and sum(self.data[self.SampleNameBox.currentText()][celllist[j]][gene].tolist()) == 0:
                    continue
                U1, p = mannwhitneyu(self.data[self.SampleNameBox.currentText()][celllist[i]][gene].tolist(), self.data[self.SampleNameBox.currentText()][celllist[j]][gene].tolist())
                corr[i,j] = p
                corr[j,i] = p
        
        sns.heatmap(corr,cmap="coolwarm_r",annot=True,xticklabels=celllist,yticklabels=celllist,)
        plt.tick_params(axis='both', which='major', labelsize=10, labelbottom = False, bottom=False, top = False, labeltop=True,labelrotation = 45)

        self.plotdict['de'] = fig
        return fig


    def drawCorr(self,subtypes):

        fig, axes = plt.subplots(1,1)
        gene = self.GeneName_Line.text()
        gene2 = self.GeneName_Line_2.text()
        if gene2 == '':
            return fig
        subdflist = []
        for i in list(self.data[self.SampleNameBox.currentText()].keys()):
                subdflist.append(self.data[self.SampleNameBox.currentText()][i].iloc[:,[0,1,2,self.data[self.SampleNameBox.currentText()][i].columns.tolist().index(gene),self.data[self.SampleNameBox.currentText()][i].columns.tolist().index(gene2)]])
        subdata = pd.concat(subdflist)
        if 'all' in subtypes:
            corr = pearsonr(subdata[gene],subdata[gene2])
            text = 'r=%s,\n p=%s' % (corr[0], corr[1])
            sns.regplot(x = gene,y =gene2,data=subdata,scatter_kws={'s':8},color='blue')
            plt.text(0.5,1.1, text, fontsize=10,horizontalalignment='center',verticalalignment='top',transform = axes.transAxes)
        else:
            subsubdata = subdata[subdata['cellType'].isin(subtypes)]
            corr = pearsonr(subsubdata[gene],subsubdata[gene2])
            text = 'r=%s,\n p=%s' % (corr[0], corr[1])
            sns.regplot(x = gene,y =gene2,data=subsubdata,scatter_kws={'s':8},color='blue')
            plt.text(0.5,1.1, text, fontsize=10,horizontalalignment='center',verticalalignment='top',transform = axes.transAxes)

        self.plotdict['cor'] = fig
        return fig


    def drawumap(self,subtypes):

        fig, axes = plt.subplots(1,1)
        gene = self.GeneName_Line.text()
        subdflist = []
        for i in list(self.data[self.SampleNameBox.currentText()].keys()):
            print(self.SampleNameBox.currentText())
            print(gene in self.data[self.SampleNameBox.currentText()][i].columns.tolist())
            subdflist.append(self.data[self.SampleNameBox.currentText()][i].iloc[:,[0,1,2,self.data[self.SampleNameBox.currentText()][i].columns.tolist().index(gene)]])
        subdata = pd.concat(subdflist)
        if 'all' in subtypes:
            zerodf = subdata[subdata[gene] == 0]
            nonzerodf = subdata[subdata[gene] != 0]
            if len(zerodf.index) != 0:
                plt.scatter(x=zerodf.iloc[:,0],y=zerodf.iloc[:,1],s=self.pointsize_box.value(),c='grey')
            if len(nonzerodf.index) != 0:
                plt.scatter(x=nonzerodf.iloc[:,0],y=nonzerodf.iloc[:,1],s=self.pointsize_box.value(),c=nonzerodf[gene], vmin=0, vmax=max(nonzerodf[gene]), cmap = 'Reds')
            plt.xlabel(subdata.columns.tolist()[0])
            plt.ylabel(subdata.columns.tolist()[1])
            plt.title(gene + ' in ' + self.SampleNameBox.currentText())
            plt.colorbar()
        else:
            subsubdata = subdata[subdata['cellType'].isin(subtypes)]
            zerodf = subsubdata[subsubdata[gene] == 0]
            nonzerodf = subsubdata[subsubdata[gene] != 0]
            if len(zerodf.index) != 0:
                plt.scatter(x=zerodf.iloc[:,0],y=zerodf.iloc[:,1],s=self.pointsize_box.value(),c='grey')
            if len(nonzerodf.index) != 0:
                plt.scatter(x=nonzerodf.iloc[:,0],y=nonzerodf.iloc[:,1],s=self.pointsize_box.value(),c=nonzerodf[gene], vmin=0, vmax=max(nonzerodf[gene]), cmap = 'Reds')
            plt.xlabel(subsubdata.columns.tolist()[0])
            plt.ylabel(subsubdata.columns.tolist()[1])
            plt.title(gene + ' in ' + self.SampleNameBox.currentText())
            plt.colorbar()

        self.plotdict['umap'] = fig
        return fig

    def drawsubumap(self,subtypes):
        fig, axes = plt.subplots(1,1)
        gene = self.GeneName_Line.text()
        subdflist = []
        for i in list(self.data[self.SampleNameBox.currentText()].keys()):
                subdflist.append(self.data[self.SampleNameBox.currentText()][i].iloc[:,[0,1,2]])
        subdata = pd.concat(subdflist)
        if 'all' in subtypes:
            groups = subdata.groupby("cellType")
            for name, group in groups:
                plt.scatter(group["UMAP_1"], group["UMAP_2"], label=name,s=self.pointsize_box.value())
            plt.legend()
            #plt.scatter(x=subdata.iloc[:,0],y=subdata.iloc[:,1],s=self.pointsize_box.value(),c=subdata.iloc[:,2])
            plt.xlabel(subdata.columns.tolist()[0])
            plt.ylabel(subdata.columns.tolist()[1])
            plt.title(gene + ' in ' + self.SampleNameBox.currentText())
        else:
            subsubdata = subdata[subdata['cellType'].isin(subtypes)]
            groups = subsubdata.groupby("cellType")
            for name, group in groups:
                plt.scatter(group["UMAP_1"], group["UMAP_2"], label=name,s=self.pointsize_box.value())
            plt.legend()
            #sns.lmplot(x='UMAP_1', y='UMAP_2', data=subsubdata, hue='cellType', fit_reg=False)
            #plt.scatter(x=subsubdata.iloc[:,0],y=subsubdata.iloc[:,1],s=self.pointsize_box.value(),c=subdata.iloc[:,2])
            plt.xlabel(subsubdata.columns.tolist()[0])
            plt.ylabel(subsubdata.columns.tolist()[1])
            plt.title(gene + ' in ' + self.SampleNameBox.currentText())

        self.plotdict['subtypeumap'] = fig
        return fig

    def drawViolin(self,subtypes):

        fig, axes = plt.subplots(1,1)
        gene = self.GeneName_Line.text()
        subdflist = []
        for i in list(self.data[self.SampleNameBox.currentText()].keys()):
                subdflist.append(self.data[self.SampleNameBox.currentText()][i].iloc[:,[0,1,2,self.data[self.SampleNameBox.currentText()][i].columns.tolist().index(gene)]])
        subdata = pd.concat(subdflist)
        if 'all' in subtypes:
            sns.violinplot(x="cellType", y=gene, data=subdata)
            plt.xlabel('Cell Type')
            plt.ylabel('Expression')
            plt.title(gene + ' in ' + self.SampleNameBox.currentText())
            plt.ylim(0)
        else:
            subdata = subdata[subdata['cellType'].isin(subtypes)]
            sns.violinplot(x="cellType", y=gene, data=subdata)
            plt.xlabel('Cell Type')
            plt.ylabel('Expression')
            plt.title(gene + ' in ' + self.SampleNameBox.currentText())
            plt.ylim(0)

        self.plotdict['violin'] = fig
        return fig

    def drawBar(self,subtypes):
        fig, axes = plt.subplots(1,1)
        gene = self.GeneName_Line.text()
        subdflist = []
        for i in list(self.data[self.SampleNameBox.currentText()].keys()):
                subdflist.append(self.data[self.SampleNameBox.currentText()][i].iloc[:,[0,1,2,self.data[self.SampleNameBox.currentText()][i].columns.tolist().index(gene)]])
        subdata = pd.concat(subdflist)
        if 'all' in subtypes:
            sns.barplot(x="cellType", y=gene, data=subdata, ci=60)
            plt.xlabel('Cell Type')
            plt.ylabel('Expression')
            plt.title(gene + ' in ' + self.SampleNameBox.currentText())
            plt.ylim(0)
        else:
            subdata = subdata[subdata['cellType'].isin(subtypes)]
            sns.barplot(x="cellType", y=gene, data=subdata, ci=60)
            plt.xlabel('Cell Type')
            plt.ylabel('Expression')
            plt.title(gene + ' in ' + self.SampleNameBox.currentText())
            plt.ylim(0)

        self.plotdict['bar'] = fig
        return fig


    def drawCGGASurvival(self):
        fig, axes = plt.subplots(1,1)
        gene = self.GeneName_Line.text()
        print(self.cggadata.columns.tolist())
        if gene not in self.cggadata.columns.tolist():
            return fig
        SurvdfU = self.cggadata.nlargest(int(self.cggadata.shape[0]/3),gene)
        SurvdfD = self.cggadata.nsmallest(int(self.cggadata.shape[0]/3),gene)

        TU = SurvdfU['T']
        EU = SurvdfU['E']
        TD = SurvdfD['T']
        ED = SurvdfD['E']
        results = logrank_test(TU,TD,EU,ED)
        #results.print_summary()
        kmf = KaplanMeierFitter()
        kmf.fit(TD,ED)
        ax = kmf.plot(label = 'Low')
        plt.hlines(y = 0.5, xmin = 0, xmax = kmf._median,linestyle = 'dotted',color='black')
        plt.vlines(x = kmf._median, ymax = 0.5, ymin = 0, linestyle = 'dotted',color='black')
        kmf.fit(TU,EU)
        ax = kmf.plot(label = 'High')
        ax.text(0.5,1.1, ('P-value: '+str(results.p_value)), fontsize=10,horizontalalignment='center',verticalalignment='top',transform = ax.transAxes)
        plt.hlines(y = 0.5, xmin = 0, xmax = kmf._median,linestyle = 'dotted',color='black')
        plt.vlines(x = kmf._median, ymax = 0.5, ymin = 0, linestyle = 'dotted',color='black')
        plt.xlim(left=0)
        plt.ylim(top=1,bottom=0)
        #plt.savefig('GBM_CGGA_NONMUT_'+i+'_Survival.pdf',format = 'pdf')
        self.plotdict['CGGASurvival'] = fig
        return fig

    def drawTCGASurvival(self):
        fig, axes = plt.subplots(1,1)
        gene = self.GeneName_Line.text()
        SurvdfU = self.tcgadata.nlargest(int(self.tcgadata.shape[0]/3),gene)
        SurvdfD = self.tcgadata.nsmallest(int(self.tcgadata.shape[0]/3),gene)

        TU = SurvdfU['T']
        EU = SurvdfU['E']
        TD = SurvdfD['T']
        ED = SurvdfD['E']
        results = logrank_test(TU,TD,EU,ED)
        #results.print_summary()
        kmf = KaplanMeierFitter()
        kmf.fit(TD,ED)
        ax = kmf.plot(label = 'Low')
        plt.hlines(y = 0.5, xmin = 0, xmax = kmf._median,linestyle = 'dotted',color='black')
        plt.vlines(x = kmf._median, ymax = 0.5, ymin = 0, linestyle = 'dotted',color='black')
        kmf.fit(TU,EU)
        ax = kmf.plot(label = 'High')
        ax.text(0.5,1.1, ('P-value: '+str(results.p_value)), fontsize=10,horizontalalignment='center',verticalalignment='top',transform = ax.transAxes)
        plt.hlines(y = 0.5, xmin = 0, xmax = kmf._median,linestyle = 'dotted',color='black')
        plt.vlines(x = kmf._median, ymax = 0.5, ymin = 0, linestyle = 'dotted',color='black')
        plt.xlim(left=0)
        plt.ylim(top=1,bottom=0)
        #plt.savefig('GBM_CGGA_NONMUT_'+i+'_Survival.pdf',format = 'pdf')
        self.plotdict['TCGASurvival'] = fig
        return fig

    # visualize
    def violin_plot(self):
        fig, axes = plt.subplots(nrows=len(self.genelist), ncols=1, figsize=(5, 4 * len(self.genelist)),
                                 constrained_layout=True)
        if len(self.genelist) == 1:
            axes = [axes]

        if self.filter_1.currentText() != '无':
            groupby = self.filter_1.currentText()
        else:
            groupby = None
        for i in range(len(self.genelist)):
            sc.pl.violin(self.data, keys=self.genelist[i], groupby=groupby, show=False, ax=axes[i])
        self.plotdict['violin'] = fig
        return fig

    def bubble_plot(self):
        fig, axes = plt.subplots()
        sc.pl.dotplot(self.data, var_names=self.genelist, groupby=self.filter_1.currentText(), show=False, ax=axes)
        self.plotdict['bubble'] = fig
        return fig

    def kaplan_meier(self):
        "Kaplan-meier model"
        figsize = (5, 4 * len(self.genelist))
        fig, axes = plt.subplots(nrows=len(self.genelist), ncols=1, figsize=figsize, constrained_layout=True)
        if len(self.genelist) == 1:
            axes = [axes]
        groupby = self.filter_1.currentText()
        group = self.comboBox_filter2.currentText()
        if group == '无':
            grouped_data = self.data.copy()
            grouped_clinical = grouped_data.obs
            title = lambda x: "Gene: %s" % (gene)
        else:
            grouped_data = self.data[self.data.obs[groupby] == group, :].copy()
            grouped_clinical = grouped_data.obs
            title = lambda x: "Group:%s Gene:%s" % (group, x)
        sample_size = int(self.horizontalSlider_boundaries.value() * grouped_clinical.shape[0] / 100)
        km_statistical_parameter = pd.DataFrame(columns=['high expression median', 'low expression median', 'p-value'])

        for i in range(len(self.genelist)):
            ax = axes[i]
            gene = self.genelist[i]
            ax.set_title(title(gene))
            # high expression
            samples_high = grouped_data.to_df()[gene].nlargest(sample_size).index
            T_high = grouped_clinical['time'].loc[samples_high]
            E_high = grouped_clinical['status'].loc[samples_high]
            kmf.fit(T_high, E_high, label='high expression')
            high_median = kmf._median
            kmf.plot(ax=ax)
            # low expression
            samples_low = grouped_data.to_df()[gene].nsmallest(sample_size).index
            T_low = grouped_clinical['time'].loc[samples_low]
            E_low = grouped_clinical['status'].loc[samples_low]
            kmf.fit(T_low, E_low, label='low expression')
            low_median = kmf._median
            kmf.plot(ax=ax)
            statistical_parameter = logrank_test(T_high, T_low, E_high, E_low, alpha=.99)
            ax.text(s='p-value:%s' % round(statistical_parameter.p_value, 4), x=grouped_clinical['time'].max() * 0.8,
                    y=0.82)
            km_statistical_parameter.loc[gene] = [high_median, low_median, statistical_parameter.p_value]
        self.statdict['km'] = km_statistical_parameter
        return fig

    def cox_regression(self):
        figsize = (6, max(4, len(self.genelist) * 0.4) + 0.5)
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
        cph = CoxPHFitter()
        cox_df = pd.concat([self.data.obs[['time', 'status']], self.data.to_df()[self.genelist]], axis=1, sort=False)
        cox_df.fillna(0, inplace=True)
        cph.fit(cox_df, 'time', 'status')
        cph.plot(ax=ax)
        return fig

    def single_cell_plot(self):
        adata = self.sc_result
        genes = ['louvain']
        genes.extend([gene for gene in self.genelist if gene in adata.var_names])
        fig, axes = plt.subplots(nrows=len(genes), ncols=1, constrained_layout=True,
                                 figsize=(5.5, 4 * len(genes)))
        if len(genes) == 1:
            axes = [axes]
        for i in range(len(genes)):
            sc.pl.umap(adata, color=genes[i], show=False, ax=axes[i])
        # self.backfeed_update('%s not include in single result.' % (
        #     ','.join([gene for gene in self.genelist if gene not in adata.var_names])))
        return fig


    def zipdir(self,path,ziph):
        for root, dirs, files in os.walk(path):
            for file in files:
                ziph.write(os.path.join(root, file))
    # output
    def save_fig(self):
        """Preview plot
        1. violin plot
        2. bubble plot
        3. survival curve
        4. cox regrssion
        5. single cell series(high expression gene,)"""

        filename = self.SampleNameBox.currentText() + '_' + self.GeneName_Line.text()

        subtypes = []
        for line in self.lines:
            if line[0].isChecked():
                subtypes.append(line[1])
                filename = filename + '_' + line[1]

        filename = QFileDialog.getSaveFileName(self.centralwidget, "保存文件",filename+'.pdf')[0]
        #filepath = folder +'/'+filename
        # os.mkdir(filepath)
        print('Saving file to:' + filename)

        # pdf = PdfPages(filepath+'/')
        # canvas = self.violin_plot()
        # pdf.savefig(canvas.figure, bbox_inches='tight')
        if filename != '':
            print('Start Saving Data...')
            pp = PdfPages(filename)
            fig = self.plotdict['umap'] 
            pp.savefig(fig)
            fig = self.plotdict['subtypeumap'] 
            pp.savefig(fig)
            fig = self.plotdict['violin']  
            pp.savefig(fig)
            fig = self.plotdict['bar'] 
            pp.savefig(fig)
            fig = self.plotdict['TCGASurvival'] 
            pp.savefig(fig)
            fig = self.plotdict['CGGASurvival'] 
            pp.savefig(fig)
            if 'cor' in list(self.plotdict.keys()):
                fig = self.plotdict['cor'] 
                pp.savefig(fig)
            fig = self.plotdict['de'] 
            pp.savefig(fig)
            pp.close()
            print('Finished Saving Data!')



        # zipf = zipfile.ZipFile(filepath+'.zip', 'w', zipfile.ZIP_DEFLATED)
        # zipdir('tmp/', zipf)
        # zipf.close()

        # pdf = PdfPages(filepath)
        # if self.UMAPbutton.isChecked():
        #     canvas = self.drawumap(subtypes)
        # elif self.Violinbutton.isChecked():
        #     canvas = self.drawViolin(subtypes)
        # pdf.savefig(canvas.figure, bbox_inches='tight')
        # pdf.close()
        # self.backfeed_update('Pdf saved.')

    def save_running_report(self):
        file_path = QFileDialog.getSaveFileName(self.centralwidget, "保存文件", ".txt", "text files (*.txt)")
        with open(file_path[0]) as f:
            f.write(self.running_report)

    def export_result(self):
        if self.tabWidget.currentWidget() == self.tab_survival and self.radioButton_km.isChecked():
            file_path = list(QFileDialog.getSaveFileName(self.centralwidget, "保存文件", ".csv", "text files (*.csv)"))[0]
            self.statdict['km'].write_csv(file_path)
        elif self.tabWidget.currentWidget() == self.tab_survival and self.radioButton_cox.isChecked():
            file_path = list(QFileDialog.getSaveFileName(self.centralwidget, "保存文件", ".csv", "text files (*.txt)"))[0]
            with open(file_path) as f:
                stdout_true = sys.stdout
                sys.stdout = f
                cph = CoxPHFitter()
                cox_df = pd.concat([self.data.obs[['time', 'status']], self.data.to_df()[self.genelist]], axis=1,
                                   sort=False)
                cox_df.fillna(0, inplace=True)
                cph.fit(cox_df, 'time', 'status')
                cph.print_summary()
                sys.stdout = stdout_true
        elif self.tabWidget.currentWidget() == self.tab_sc:
            file_path = \
                list(QFileDialog.getSaveFileName(self.centralwidget, "保存文件", ".h5ad", "Anndata H5 file (*.h5ad)"))[0]
            self.sc_result.write_h5ad(file_path)
        else:
            self.backfeed_update('当前分析标签页没有分析结果.')

    def single_cell_start(self):
        self.backfeed_update('Single cell analysis start.')
        from platform_childthread import SingleCell
        parameters = [self.lineEdit_min_genes_cell.text(),
                      self.lineEdit_min_counts_cell.text(),
                      self.lineEdit_max_counts_cell.text(),
                      self.lineEdit_max_mito_cell.text(),
                      self.lineEdit_min_cells_gene.text(),
                      self.lineEdit_min_counts_gene.text(),
                      self.lineEdit_pcs.text()]
        parameters = tuple([int(i) for i in parameters])
        print(parameters)
        sc_thread = SingleCell(self.data, parameters)
        sc_thread.signal.connect(self.single_cell_end)
        self.threadpool['SingleCell'] = sc_thread
        sc_thread.start()

    def single_cell_end(self):
        sc_result = self.threadpool['SingleCell'].sc_result
        if type(sc_result) == str:
            self.backfeed_update(sc_result)
        else:
            self.sc_result = self.threadpool['SingleCell'].sc_result
            if sc_result.obsm['X_pca'].shape[1] < int(self.lineEdit_pcs.text()):
                self.backfeed_update('X_pca does not have enough PCs. The value "Size of PCs" have reset to %s' %
                                     sc_result.obsm['X_pca'].shape[1])
            sc.pl.violin(sc_result, ['total_counts', 'n_genes_by_counts', 'mito_percent'], multi_panel=True)
            sc.pl.scatter(sc_result, x='total_counts', y='mito_percent')
            sc.pl.scatter(sc_result, x='total_counts', y='n_genes_by_counts')
            sc.pl.highest_expr_genes(sc_result)
            sc.pl.pca_variance_ratio(sc_result)
            plt.close()
            self.backfeed_update('Single cell analysis finished.')
            self.backfeed_update(sc_result.__str__())


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    app.setStyle('Fusion')
    MainWindow = QtWidgets.QMainWindow()
    ui = myGUI()
    ui.setupUi(MainWindow)
    ui.setup_signal()
    ui.set_default()
    MainWindow.show()
    sys.exit(app.exec_())
