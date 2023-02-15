#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
    Workflow for single cell RNA sequencing analysis.
    
@author: earezza
"""
import itertools
import scanpy as sc
import seaborn as sns
import pandas as pd
import numpy as np
import anndata as ad
import palantir
import scvelo as scv
import matplotlib.pyplot as plt
'''
from bokeh.plotting import figure, output_file, show, curdoc
from bokeh.models import GlyphRenderer, Circle, ColumnDataSource, CDSView, IndexFilter, BooleanFilter, GroupFilter, Button, Panel, Tabs, Legend, LegendItem, ColorBar, CustomJS, Select, MultiChoice, TextInput, CategoricalColorMapper
from bokeh.layouts import gridplot, column, row
from bokeh.palettes import Set1_6, Category20_20
from bokeh.transform import factor_cmap, factor_mark, linear_cmap, log_cmap
from bokeh.io import output_notebook, curdoc
from bokeh.themes import Theme
'''
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
sc.set_figure_params(scanpy=True, dpi=300, dpi_save=300, 
                     frameon=True, vector_friendly=True, 
                     fontsize=14, figsize=None, 
                     color_map=None, format='pdf', 
                     facecolor=None, transparent=False, 
                     ipython_format='png2x')
#plt.rcParams['figure.dpi'] = 300
#plt.rcParams['savefig.dpi'] = 300
#sns.set(rc={"figure.dpi":300, 'savefig.dpi':300, 'axes.facecolor':'white', 'figure.facecolor':'white'})

# organism as “hsapiens”, “mmusculus”, “drerio”, etc.
# host as “grch37.ensembl.org” or valid BioMart URLs
ensembl_annot = sc.queries.biomart_annotations('mmusculus', ["ensembl_gene_id", 'external_gene_name', 'hgnc_symbol', "start_position", "end_position", "chromosome_name"], host='www.ensembl.org')

# Taken from R package Seurat cc.genes
cell_cycle_genes = {'S': ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
                          "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1",
                          "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2",
                          "CDC45", "CDC6" , "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1",
                          "CHAF1B", "BRIP1", "E2F8"
                          ],
                    'G2M': ["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B",
                            "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1",
                            "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C",
                            "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1",
                            "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"
                            ]
                    }
cell_cycle_genes_update2019 = {'S': ["MCM5", "PCNA", "TYMS", "FEN1", "MCM7", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
                          "CDCA7", "DTL", "PRIM1", "UHRF1", "CENPU", "HELLS", "RFC2", "POLR1B", "NASP", "RAD51AP1",
                          "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2",
                          "CDC45", "CDC6" , "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1",
                          "CHAF1B", "MRPL36", "E2F8"
                          ],
                    'G2M': ["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B",
                            "MKI67", "TMPO", "CENPF", "TACC3", "PIMREG", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1",
                            "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "JPT1", "CDC20", "TTK", "CDC25C",
                            "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1",
                            "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"
                            ]
                    }
# Format according to data's gene naming
for i in range(0, len(cell_cycle_genes['S'])):
    cell_cycle_genes['S'][i] = cell_cycle_genes['S'][i].lower().capitalize()
for i in range(0, len(cell_cycle_genes['G2M'])):
    cell_cycle_genes['G2M'][i] = cell_cycle_genes['G2M'][i].lower().capitalize()
for i in range(0, len(cell_cycle_genes_update2019['S'])):
    cell_cycle_genes_update2019['S'][i] = cell_cycle_genes_update2019['S'][i].lower().capitalize()
for i in range(0, len(cell_cycle_genes_update2019['G2M'])):
    cell_cycle_genes_update2019['G2M'][i] = cell_cycle_genes_update2019['G2M'][i].lower().capitalize()
    
# Change according to data analysis...e.g. 'Cell Type': ['list', 'of', 'gene', 'markers']
markers = {'MSC': ['Lepr', 'Cxcl12', 'Adipoq', 'Kitl', 'Angpt1', 'Nt5e', 'Vcam1', 'Eng', 'Grem1'],
           'OLC': ['Bglap'],
           'Chondrocyte': ['Col2a1', 'Acan'],
           'Fibroblast': ['S100a4'],
           'BMEC': ['Cdh5'],
           'Pericyte': ['Acta2']
          }
genes = list(itertools.chain(*list(markers.values())))

ko = sc.read_10x_mtx('filtered_feature_bc_matrix_Ko/', cache=True)
wt = sc.read_10x_mtx('filtered_feature_bc_matrix_Wt/', cache=True)                              
ko.var_names_make_unique()  
wt.var_names_make_unique() 

ko.obs['sample'] = "KO"
wt.obs['sample'] = "WT"
data = ko.concatenate(wt)
del ko, wt  # delete no longer used variables to save memory  

data
data.obs
data.var

data.var['mt'] = data.var_names.str.startswith('MT-')
data.var['ribo'] = data.var_names.str.startswith('ribo-')

sc.pp.calculate_qc_metrics(data, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
data

sns.jointplot(
    data=data.obs,
    x="total_counts",
    y="n_genes_by_counts",
    kind="reg",
    marginal_ticks=True,
    color='indigo',
    height=7, # changes size of plot
).fig.suptitle("%s cells"%data.obs.shape[0])

sns.scatterplot(data=data.obs, x='total_counts', y='pct_counts_mt', color='indigo')
sns.scatterplot(data=data.obs, x='total_counts', y='pct_counts_ribo', color='indigo')

sc.pl.violin(data, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True, color='indigo')

data = data[data.obs.pct_counts_mt < 5, :] # keeps all cells where < 5% of expression counts are mitochondrial genes
data = data[data.obs.pct_counts_ribo > 5, :] # keeps all cells where > 5% of expression counts are ribosomal genes

min_num_genes = data.obs['n_genes_by_counts'].quantile(q=0.01, interpolation='lower')
max_num_genes = data.obs['n_genes_by_counts'].quantile(q=0.99, interpolation='higher')
print('Minimum number of genes = %s'%min_num_genes)
print('Maxmium number of genes = %s'%max_num_genes)

sc.pp.filter_cells(data, min_genes=min_num_genes) # Filter out low gene expressions
sc.pp.filter_cells(data, max_genes=max_num_genes) # Filter out high gene expressions

sc.pp.filter_genes(data, min_cells=3) # Filter out genes with < 3 cells expressing them

data = data[:, ~data.var_names.isin(cell_cycle_genes['S'])]
data = data[:, ~data.var_names.isin(cell_cycle_genes['G2M'])]
data = data[:, ~data.var_names.isin(cell_cycle_genes_update2019['S'])]
data = data[:, ~data.var_names.isin(cell_cycle_genes_update2019['G2M'])]

sc.pl.violin(data, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True, color='indigo')

data

sc.pp.normalize_total(data, target_sum=1e4)
sc.pp.log1p(data)

# sc.pp.highly_variable_genes(data, flavor='seurat', min_mean=0.0125, max_mean=3, min_disp=0.5) # Use if single sample loaded in dataset
sc.pp.highly_variable_genes(data, flavor='seurat', min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key='sample') # Use if multiple samples loaded in dataset
sc.pl.highly_variable_genes(data)

var_genes = data.var[data.var['highly_variable_nbatches'] > 1] # highly variable among >1 batches (here, >1 is both samples, i.e. nbatches=2)
var_genes

data
data.raw = data
data.write('data_clean.h5ad')

data = sc.read('data_clean.h5ad')

data_batch_corrected = sc.AnnData(X=data.raw.X, var=data.raw.var, obs = data.obs)
data_batch_corrected.raw = data_batch_corrected

sc.pp.combat(data_batch_corrected, key='sample')

# Other batch correction methods, if used instead of ComBat, these may require additional package installation
#sc.external.pp.bbknn(data_batch_corrected, key='sample')
#sc.external.pp.harmony_integrate(data_batch_corrected, key='sample')
#sc.external.pp.mnn_correct(data_batch_corrected, key='sample')
#sc.external.pp.scanorama_integrate(data_batch_corrected, key='sample')

data_batch_corrected.raw = data_batch_corrected
sc.pp.scale(data_batch_corrected, max_value=10)

data_batch_corrected.write('data_batch_corrected.h5ad')

data = sc.read('data_batch_corrected.h5ad')
data

sc.tl.pca(data, n_comps=50, svd_solver='arpack', random_state=2022)
sc.pl.pca(data)
sc.pl.pca_variance_ratio(data, log=True)

sc.pp.neighbors(data, knn=True, n_neighbors=10, n_pcs=40, random_state=2022)

sc.tl.tsne(data, n_pcs=50, use_fast_tsne=False, random_state=2022) # May take a minute...
sc.tl.umap(data, random_state=2022) # recommended
sc.tl.draw_graph(data, random_state=2022)

sc.tl.louvain(data, resolution=1.0, random_state=2022) # recommended
sc.tl.leiden(data, resolution=1.0, random_state=2022)

data

sc.pl.tsne(data, color=["louvain"], legend_loc='on data', color_map='viridis')
sc.pl.umap(data, color=["louvain"], legend_loc='on data', color_map='viridis')
sc.pl.draw_graph(data, color=["louvain"], legend_loc='on data', color_map='viridis')

genes = ['provide', 'list', 'of', 'genes', 'to', 'see', 'here']

sc.pl.tsne(data, color=genes, legend_loc='on data', color_map='viridis')
sc.pl.umap(data, color=genes, legend_loc='on data', color_map='viridis')
sc.pl.draw_graph(data, color=genes, legend_loc='on data', color_map='viridis')

sc.pl.violin(data, genes, group='louvain')

'''
clusters = ['2', '3', '4', '5', '7', '9', '13', '15'] # Change values to desired clusters
sub_data = data[data.obs['louvain'].isin(clusters)]
# Extract matrix
adata = sc.AnnData(sub_data.X)
adata.obs = sc.get.obs_df(sub_data)
adata.var = sc.get.var_df(sub_data)
# Keep global clustering annotations
adata.obs['louvain_global'] = sub_data[sub_data.obs_names.isin(adata.obs_names)].obs['louvain']

# Save to file
adata.write_h5ad('sub_data.h5ad')
'''

'''
# Load data
#data = sc.read('data_batch_corrected.h5ad')
data = sc.read('yuefeng/write/data_batch_corrected.h5ad')

# Extract DataFrame object from matrix info
df = pd.DataFrame(data.raw.X, columns=data.var_names, index=data.obs_names)
df.insert(df.shape[1], 'UMAP1', data.obsm['X_umap'][:,0])
df.insert(df.shape[1], 'UMAP2', data.obsm['X_umap'][:,1])
df = pd.merge(df, data.obs, left_index=True, right_index=True)

# Load df subset for Bokeh example
df = pd.read_csv('obs_sample.tsv', sep='\t', index_col=0)

# Define interactive Bokeh plot to be displayed
def scviewer(doc):
    LOUVAIN_GROUPS = np.sort(df['louvain'].unique()).astype(str)
    GENES = sorted(list(df.columns[:-10].values))
    VIEW_OPTIONS = ['None', 'louvain', 'sample'] + GENES
    SAMPLES = list((set(df['sample'].unique())))
    sample_colors = factor_cmap('sample', 'Set1_6', factors=SAMPLES)
    louvain_colors = factor_cmap('louvain', 'Category20_20', factors=LOUVAIN_GROUPS)

    df['louvain'] = df['louvain'].astype(str)
    source = ColumnDataSource(data=df)

    TOOLTIPS = [
        ("louvain", "@louvain"),
        ("sample", "@sample"),
        ("total_counts", "@total_counts"),
        ("n_genes", "@n_genes")
    ]

    # Update plot colors for each cell
    def update_color(new):
        # Default
        if new == 'None':
            renders['main'].glyph.fill_color = 'darkgrey'
            renders['main'].visible = True
            for group in SAMPLES:
                renders[group].visible = False
            for group in LOUVAIN_GROUPS:
                renders[group].visible = False
            color_bar.visible = False

        # Color by sample
        elif new == 'sample':
            for group in SAMPLES:
                renders[group].visible = True
            for group in LOUVAIN_GROUPS:
                renders[group].visible = False
            renders['main'].visible = False
            color_bar.visible = False

        # Color by gene expression
        elif new in GENES:
            renders['main'].glyph.fill_color = linear_cmap(new, 'Viridis256', source.data[new].min(), source.data[new].max())
            renders['main'].visible = True
            color_bar.visible = True
            for group in SAMPLES:
                renders[group].visible = False
            for group in LOUVAIN_GROUPS:
                renders[group].visible = False

        # Color by louvain group
        elif new == 'louvain':
            for group in LOUVAIN_GROUPS:
                renders[group].visible = True
            for group in SAMPLES:
                renders[group].visible = False
            renders['main'].visible = False
            color_bar.visible = False

        else:
            renders['main'].glyph.fill_color = 'black'
            renders['main'].visible = True
            for group in SAMPLES:
                renders[group].visible = False
            for group in LOUVAIN_GROUPS:
                renders[group].visible = False
            color_bar.visible = False

    # Update plot legend based on selected view

    #def update_legend(new):
    #    if new == 'sample':
    #        sample_legend.visible = True
    #        louvain_legend.visible = False
    #    if new == 'louvain':
    #        louvain_legend.visible = True
    #        sample_legend.visible = False
    #    else:
    #        sample_legend.visible = False
    #        louvain_legend.visible = False

    # Handler for updating plot
    def update_plot(attr, old, new):
        update_color(new)
        #update_legend(new)


    # Create plot
    fig = figure(title='YL Data (n_cells=%s n_genes=%s)'%(df.shape[0], len(GENES)), 
                 toolbar_location="above",
                 tools="hover,pan,box_zoom,wheel_zoom,lasso_select,box_select,save,undo,redo,reset", 
                 tooltips=TOOLTIPS,
                 #lod_factor=10,
                 #lod_threshold=1000,
                 sizing_mode='stretch_both',
                 )
    fig.toolbar.autohide = False
    fig.xaxis.axis_label = 'UMAP1'
    fig.yaxis.axis_label = 'UMAP2'

    # Add cells onto UMAP projection
    renders = {}
    renders['main'] = fig.circle(x='UMAP1', y='UMAP2', source=source, 
                                 alpha=1.0, fill_color='darkgrey', line_color=None, 
                                 muted_alpha=0.1, muted_fill_color='gray', visible=True)
    # Add color bar for gene expression
    color_bar = ColorBar(title='Gene Expression', #title_text_align='center', title_text_baseline='middle',
                         minor_tick_line_alpha=0, major_tick_line_alpha=0,
                         color_mapper=factor_cmap('colorbar', 'Viridis256', factors=['Low', 'High'])['transform'], 
                         visible=False)

    # Add all other renderers
    # Add renderers for samples
    samples_render = []
    sample_legend_items = []
    for group, color in zip(df['sample'].unique(), Set1_6):
        cds = ColumnDataSource(data=df[df['sample'] == group])
        renders[group] = fig.circle(x='UMAP1', y='UMAP2', source=cds, 
                                    fill_alpha=0.6, fill_color=color, line_color=color, line_alpha=1.0, 
                                    muted_alpha=0.1, muted_fill_color='grey',
                                    visible=False)
        sample_legend_items.append(LegendItem(label=group, renderers=[fig.renderers[-1]]))
    sample_legend = Legend(items=sample_legend_items, title='Sample', click_policy='hide', location='top_left')

    # Add renderers for louvain groups
    louvain_render = []
    louvain_legend_items = []
    for group, color in zip(LOUVAIN_GROUPS, Category20_20):
        cds = ColumnDataSource(data=df[df['louvain'] == group])
        renders[group] = fig.circle(x='UMAP1', y='UMAP2', source=cds, 
                                    fill_alpha=1.0, fill_color=color, line_color=None, 
                                    muted_alpha=0.1, muted_fill_color='grey',
                                    visible=False)
        louvain_legend_items.append(LegendItem(label=group, renderers=[fig.renderers[-1]]))
    louvain_legend = Legend(items=louvain_legend_items, title='Louvain', click_policy='hide', location='top_left')

    del cds

    # Add legends
    fig.add_layout(sample_legend, 'left')
    fig.add_layout(louvain_legend, 'left')
    fig.add_layout(color_bar, 'left')

    # Add dropdown selection for views
    select = Select(title="View:", value='None', options=VIEW_OPTIONS, sizing_mode='fixed', width=200, height=50)
    select.on_change('value', update_plot)

    layout = column(select, fig, sizing_mode='stretch_both')

    doc.add_root(layout)
    doc.title = 'scRNAseq Plot'
    
output_notebook() # To embed/display Bokeh server plot in this notebook
show(scviewer, notebook_url="http://localhost:8888") # May take a moment to load...
'''

sc.tl.rank_genes_groups(data, 'louvain', method='wilcoxon', random_state=2022)

ranked_genes_in_clusters = pd.DataFrame(data.uns['rank_genes_groups']['names'])
result = data.uns['rank_genes_groups']
groups = result['names'].dtype.names
groups_pval = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']})

ranked_genes_in_clusters.head(n=10) # See highest top 10 expressed genes for each cluster

groups_pval.head(n=10) # See highest top 10 expressed genes for each cluster with p-values & details

groups_pval.to_csv('ranked_genes.csv', header=True, index=True, sep=',')


sc.external.tl.palantir(data,  n_components=5, knn=30)
dm_res = {'T': data.obsp['palantir_diff_op'].T, 
          'EigenVectors': pd.DataFrame(data.obsm['X_palantir_diff_comp'], index=data.obs_names), 
          'EigenValues': pd.Series(data.uns['palantir_EigenValues']), 
          'kernel': data.obsp['palantir_diff_op']}
ms_data = pd.DataFrame(data.obsm['X_palantir_multiscale'], index=data.obs_names)

data

sc.tl.tsne(data, n_pcs=2, use_rep='X_palantir_multiscale', perplexity=150, random_state=2022)

sc.pl.tsne(data, layer='palantir_imp', color=["louvain"], legend_loc='on data', color_map='viridis')
sc.pl.umap(data, layer='palantir_imp', color=["louvain"], legend_loc='on data', color_map='viridis')
sc.pl.draw_graph(data, layer='palantir_imp', color=["louvain"], legend_loc='on data', color_map='viridis')

genes = ['provide', 'list', 'of', 'genes', 'to', 'see', 'here']

sc.pl.tsne(data, layer='palantir_imp', color=genes, legend_loc='on data', color_map='viridis')
sc.pl.umap(data, layer='palantir_imp', color=genes, legend_loc='on data', color_map='viridis')
sc.pl.draw_graph(data, layer='palantir_imp', color=genes, legend_loc='on data', color_map='viridis')

tsne = pd.DataFrame(data.obsm['X_tsne'], index=data.obs_names)
tsne.rename(columns={0: 'x', 1: 'y'}, inplace=True)

umap = pd.DataFrame(data.obsm['X_umap'], index=data.obs_names)
umap.rename(columns={0: 'x', 1: 'y'}, inplace=True)

fa = pd.DataFrame(data.obsm['X_draw_graph_fa'], index=data.obs_names)
fa.rename(columns={0: 'x', 1: 'y'}, inplace=True)

palantir.plot.plot_diffusion_components(tsne, dm_res)
palantir.plot.plot_diffusion_components(umap, dm_res)
palantir.plot.plot_diffusion_components(fa, dm_res)

early_cell = pd.DataFrame(data[:, ['Sox9']].layers['palantir_imp'], index=data.obs_names, columns=['Sox9']).idxmax().values[0] # for max expression
#early_cell = pd.DataFrame(data[:, ['Sox9']].layers['palantir_imp'], index=data.obs_names, columns=['Sox9']).idxmin().values[0] # for min expression

#terminal_cells = pd.DataFrame(data[:, ['Sox9']].layers['palantir_imp'], index=data.obs_names, columns=['Sox9']).idxmax().values # for max expression
#terminal_cell = pd.DataFrame(data[:, ['Sox9']].layers['palantir_imp'], index=data.obs_names, columns=['Sox9']).idxmin().values # for min expression

num_waypoints = int(0.05*data.n_obs)

pr_res = sc.external.tl.palantir_results(data, early_cell, ms_data='X_palantir_multiscale', num_waypoints=num_waypoints, use_early_cell_as_start=True)#, terminal_states=terminal_states.index)

palantir.plot.plot_palantir_results(pr_res, tsne)
palantir.plot.plot_palantir_results(pr_res, umap)
palantir.plot.plot_palantir_results(pr_res, fa)

cells = pr_res.branch_probs.columns.values.tolist()
cells.append(early_cell)
palantir.plot.plot_terminal_state_probs(pr_res, cells)

palantir.plot.highlight_cells_on_tsne(tsne, cells, start_cell=early_cell)
palantir.plot.highlight_cells_on_tsne(umap, cells, start_cell=early_cell)
palantir.plot.highlight_cells_on_tsne(fa, cells, start_cell=early_cell)

genes = ['provide', 'list', 'of', 'genes', 'to', 'see', 'here']

imp_df = pd.DataFrame(data[:, genes].layers['palantir_imp'], index=data.obs_names, columns=genes) # Store imputed data into DataFrame
gene_trends = palantir.presults.compute_gene_trends(pr_res, imp_df.loc[:, genes])

palantir.plot.plot_gene_trends(gene_trends)
palantir.plot.plot_gene_trend_heatmaps(gene_trends)

imp_df_first1000 = pd.DataFrame(data[:, data.var_names.values[0:1000]].layers['palantir_imp'],
                     index=data.obs_names, columns=data.var_names.values[0:1000]) # Store imputed data into DataFrame
gene_trends_first1000 = palantir.presults.compute_gene_trends(pr_res, imp_df_first1000.iloc[:, 0:1000], )

for c in range(0, len(cells)-1):
    trends = gene_trends_first1000[cells[c]]['trends']
    gene_clusters = palantir.presults.cluster_gene_trends(trends)
    palantir.plot.plot_gene_trend_clusters(trends, gene_clusters)
    
    

