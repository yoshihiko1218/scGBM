from PyQt5.Qt import QThread, pyqtSignal, QWidget, QMessageBox,QApplication,QPushButton,QVBoxLayout
from PyQt5 import QtWidgets
import anndata as ad
import sys
import logging
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import pickle


class Loading(QThread):
    signal = pyqtSignal()

    def __init__(self, filepath, tcgapath, cggapath):
        super().__init__()
        self.filepath = filepath
        self.tcgapath = tcgapath
        self.cggapath = cggapath

    def run(self):
        print('Start Loading Data....')
        with open(self.filepath, 'rb') as f:
            data = pickle.load(f)
        self.tcgadata = pd.read_csv(self.tcgapath,sep='\t',index_col=0)
        self.cggadata = pd.read_csv(self.cggapath,sep='\t',index_col=0)
        #data = CustomUnpickler(open(self.filepath, 'rb')).load()
        self.data = data
        print('Finished Loading Data!')
        self.signal.emit()

class SingleCell(QThread):
    signal = pyqtSignal()

    def __init__(self, adata, parameters):
        super().__init__()
        self.data = adata.copy()
        self.parameters = parameters

    def run(self):
        print('sc start')
        adata = self.data.copy()
        shape_list = []
        min_genes, min_counts_cell, max_counts_cell, max_mito, min_cells, min_counts_gene, n_pcs = self.parameters
        qc_df = sc.pp.calculate_qc_metrics(adata)
        adata.obs = pd.concat([adata.obs, qc_df[0]], axis=1)
        adata.var = pd.concat([adata.var, qc_df[1]], axis=1)
        mito_gene = [gene for gene in adata.var_names if gene.startswith('mt-')]
        adata.obs['mito_percent'] = adata.to_df()[mito_gene].sum(1) / adata.obs['total_counts']

        # Data Clean
        adata = adata[adata.obs['n_genes_by_counts'] >= min_genes, :]
        shape_list.append(adata.shape[0] * adata.shape[1])
        adata = adata[adata.obs['total_counts'] >= min_counts_cell, :]
        shape_list.append(adata.shape[0] * adata.shape[1])
        adata = adata[adata.obs['total_counts'] <= max_counts_cell, :]
        shape_list.append(adata.shape[0] * adata.shape[1])
        adata = adata[adata.obs['mito_percent'] < max_mito, :]
        shape_list.append(adata.shape[0] * adata.shape[1])
        adata = adata[:, adata.var['n_cells_by_counts'] >= min_cells]
        shape_list.append(adata.shape[0] * adata.shape[1])
        adata = adata[:, adata.var['total_counts'] >= min_counts_gene]
        shape_list.append(adata.shape[0] * adata.shape[1])
        print('sc filtered')
        if adata.shape[0] * adata.shape[1] > 0:
            # Normalization
            sc.pp.log1p(adata)
            sc.pp.normalize_total(adata)
            sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
            print('sc normalized')
            # regress out
            sc.pp.regress_out(adata, ['total_counts', 'mito_percent'])
            sc.pp.scale(adata, max_value=10)
            print('sc scaled')

            # reduce dimensional
            sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True)
            if adata.obsm['X_pca'].shape[1] < n_pcs:
                n_pcs = adata.obsm['X_pca'].shape[1]
            sc.pp.neighbors(adata, n_neighbors=10, n_pcs=n_pcs)
            sc.tl.umap(adata, random_state=0)
            print('sc reduced')
            # Cluster
            sc.tl.louvain(adata)
            self.sc_result = adata
            print('sc finished')
        else:
            parameter_names = ['Minmum genes per cell', 'Minmum counts per cell', 'Maximum counts per cell',
                               'Maximum mitomchondria(%)', 'Minmum cells per gene', 'Minmum counts per gene']

            self.sc_result = 'No data included in matrix after filter. Reset the "%s" value and retry.' % \
                             parameter_names[len([i for i in shape_list if i > 0])]
        self.signal.emit()

