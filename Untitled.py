#!/usr/bin/env python
# coding: utf-8

# In[97]:


get_ipython().system('pip install nbconvert')
get_ipython().system('pip install nbconvert[webpdf]')
get_ipython().system('jupyter nbconvert --to webpdf --allow-chromium-download untitled.ipynb')


# In[2]:


get_ipython().system('pip install scanpy')


# In[3]:


import scanpy as sc


# In[4]:


adata = sc.read_csv('/Users/josephanthony/Desktop/scrnaseqproject/GSM5226574_C51ctr_raw_counts.csv').T
adata


# In[5]:


adata.X.shape


# In[6]:


#Preprocessing-- Doublet Removal


# In[7]:


get_ipython().system('pip install scvi_tools')
import scvi


# In[8]:


get_ipython().system('pip3 install --user scikit-misc')
import skmisc.loess


# In[9]:


sc.pp.filter_genes(adata, min_cells = 10)


# In[10]:


adata


# In[11]:


sc.pp.highly_variable_genes(adata, n_top_genes = 2000, subset = True, flavor = 'seurat_v3')


# In[12]:


get_ipython().system('pip install --user scikit-misc')


# In[12]:


get_ipython().system('pip install scikit-misc')


# In[13]:


adata.var


# In[14]:


adata.var[adata.var.index.str.startswith('MT-')]


# In[17]:


adata.var['mt'] = adata.var.index.str.startswith('MT-')


# In[18]:


adata.var


# In[19]:


#Now for ribosomal genes


# In[20]:


import pandas as pd


# In[21]:


ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"


# In[22]:


ribo_genes = pd.read_table(ribo_url, skiprows=2, header=None)


# In[23]:


ribo_genes


# In[24]:


adata.var['ribo'] = adata.var_names.isin(ribo_genes[0].values)


# In[27]:


adata.var


# In[28]:


sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=True)


# In[30]:


adata.var.sort_values('n_cells_by_counts')


# In[34]:


adata.obs.sort_values('total_counts')


# In[32]:


sc.pp.filter_cells(adata, min_genes=200)


# In[33]:


adata.var


# In[35]:


sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'], jitter=0.4, multi_panel=True)


# In[36]:


import numpy as np


# In[37]:


upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)


# In[38]:


upper_lim


# In[39]:


adata.obs


# In[40]:


adata = adata[adata.obs.pct_counts_mt < 20]


# In[41]:


adata = adata[adata.obs.pct_counts_ribo < 2]


# In[42]:


adata


# In[43]:


#Normalization - so we can compare cells and compare genes


# In[44]:


adata.X.sum(axis=1)


# In[45]:


sc.pp.normalize_total(adata,target_sum=1e4) #normalize every cell to 10k UMI


# In[46]:


adata.X.sum(axis=1)


# In[47]:


sc.pp.log1p(adata) #change to log counts


# In[48]:


adata.X.sum(axis=1)


# In[49]:


adata.raw = adata


# In[50]:


#Clustering


# In[52]:


sc.pp.highly_variable_genes(adata, n_top_genes = 1944)


# In[54]:


adata.var


# In[55]:


sc.pl.highly_variable_genes(adata)


# In[56]:


adata = adata[:, adata.var.highly_variable]


# In[57]:


adata


# In[58]:


sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt', 'pct_counts_ribo'])


# In[59]:


sc.pp.scale(adata, max_value=10)


# In[60]:


sc.tl.pca(adata, svd_solver='arpack')


# In[61]:


sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)


# In[62]:


sc.pp.neighbors(adata, n_pcs=30)


# In[63]:


adata


# In[64]:


adata.obsp['connectivities']


# In[65]:


sc.tl.umap(adata)


# In[66]:


sc.pl.umap(adata)


# In[67]:


get_ipython().system('pip install leidenalg')


# In[68]:


sc.tl.leiden(adata, resolution=0.5)


# In[69]:


adata.obs


# In[70]:


sc.pl.umap(adata, color=['leiden'])


# In[73]:


sc.tl.rank_genes_groups(adata, "leiden")


# In[74]:


sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)


# In[75]:


markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
markers


# In[84]:


adata.obs.head()


# In[85]:


adata


# In[86]:


adata.layers['counts'] = adata.X.copy()


# In[ ]:




