# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 14:49:32 2021

@author: Daniel Keitley
"""

from . import plot_utils

import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import plotly.graph_objects as go





def createAnnotationDirs(base_path, clusters):
    """
    
    e.g.
    createAnnotationDirs(base_path = "../figs/celltype_annotation/annotation_pipeline/r_mesoderm/", 
                     clusters = ["leiden_res8","mesoderm_leiden_res2","mesoderm_leiden_res5"])
                     
    """
    
    # Marker expression plots
    Path(base_path + "marker_expression").mkdir(parents=True, exist_ok=True)
    
    # Cluster fraction plots
    Path(base_path + "cluster_fractions").mkdir(parents=True, exist_ok=True)
    
    # Prediction fraction plots
    Path(base_path + "prediction_fractions").mkdir(parents=True, exist_ok=True)
    
    for cluster in clusters:
        # Create directory for each clustering
        Path(base_path + "cluster_fractions/" + cluster).mkdir(parents=True, exist_ok=True)
        Path(base_path + "prediction_fractions/" + cluster).mkdir(parents=True, exist_ok=True)
        
    
    

def computeObsFraction(adata, obsA="leiden_res8",obsB="singler", obsB_ignore_thresh = None ):
    """
    Computes the fraction of observations B within groups of observation A.
    
    """
    
    df = pd.DataFrame({obsA:adata.obs[obsA], obsB:adata.obs[obsB]},index=adata.obs.index)
    df[obsA + "_ncells"] = df.groupby([obsA]).transform('count')

    df[obsB + "_count"] = 0

    counts = df.groupby([obsA,obsB,obsA + "_ncells"]).count()
    counts[obsB + "_frac"] = counts[obsB + "_count"].groupby(obsA).transform(lambda x: x/x.sum())
    counts.sort_values([obsA,obsB,obsB + "_frac"],ascending=False)
    
    counts= counts.reset_index()
    counts = counts.dropna()
    
     # Group observations B with very few cells as 'Other'
    if(obsB_ignore_thresh is not None) :
        counts[obsB] = counts[obsB].astype(str)
        counts.loc[counts[obsB + "_frac"] <= obsB_ignore_thresh,obsB] = "Other"      
        
    
    
    return(counts)




def plotObsFraction(data, obA, obsA, obsB, obsA_ignore_thresh = None, 
                    obsB_ignore_thresh = None, export_dir=None, 
                    obsB_colours=None,ax=None,
                   figsize=(8,8)):
    """
    obA: A specific value of obsA
    
    """
    
    # Get obs fractions if anndata object passed in
    if(isinstance(data,sc.AnnData)):
        data = computeObsFraction(data, obsA = obsA, obsB = obsB, obsB_ignore_thresh = obsB_ignore_thresh)
    
    # Don't plot fractions for observations with very few cells
    if(obsA_ignore_thresh is not None):
        data = data.loc[data[obsA + "_ncells"] >= obsA_ignore_thresh]
        
    if(isinstance(obA,list)):
            nobs = len(obA)
            df_plot = data[data[obsA].isin(obA)]
    else:
        df_plot = data[data[obsA]==obA]
    
    if(isinstance(obsB_colours,dict)):
                obsB_colours = [obsB_colours[i] for i in df_plot[obsB]]
       
    if(ax is None):
        if(isinstance(obA,list)):
            
            fig, ax = plt.subplots(1,nobs,figsize=figsize)
            
            for i in range(nobs):
                df_subplot = df_plot[df_plot[obsA]==obA[i]]
                ax[i].barh(df_subplot[obsB], df_subplot[obsB + "_frac"], color=obsB_colours)
                ax[i].set_title(obsA + ": " +  obA[i] + " (" + str(int(df_subplot[obsA + "_ncells"].iloc[0]))+ " cells)")
                ax[i].set_xlabel("Fraction of cells")
                ax[i].set_ylabel(obsB)
                plt.tight_layout()
            return(ax)
        else:
            fig, ax = plt.subplots(figsize=figsize)


    ax.barh(df_plot[obsB], df_plot[obsB + "_frac"], color=obsB_colours)
    ax.set_title(obsA + ": " +  obA + " (" + str(int(df_plot[obsA + "_ncells"].iloc[0]))+ " cells)")
    ax.set_xlabel("Fraction of cells")
    ax.set_ylabel(obsB)
            
    plt.tight_layout()
    
    return(ax)
        
    


def plotMultiObsFraction(data, obsA, obsB, obsA_ignore_thresh = None, obsB_ignore_thresh = None, export_dir=None, obsB_colours=None):
    
    """
    
    e.g. 
    plotMultiObsFraction(r_mesoderm, obsA = "leiden_res8", obsB = "singler", 
                     obsA_ignore_thresh = 50 , obsB_ignore_thresh = 0.02,
                     #colours = scrabbit.getCelltypeColours(),
                     export_dir = "../figs/celltype_annotation/annotation_pipeline/r_mesoderm/cluster_fractions/leiden_res8/")
                     
    """
    # Get obs fractions if anndata object passed in
    if(isinstance(data,sc.AnnData)):
        data = computeObsFraction(data, obsA = obsA, obsB = obsB, obsB_ignore_thresh = obsB_ignore_thresh)
    
    # Don't plot fractions for observations with very few cells
    if(obsA_ignore_thresh is not None):
        data = data.loc[data[obsA + "_ncells"] >= obsA_ignore_thresh]
    
    
    # if(export_dir is None):
        # Plot in a grid
    
    
    if(export_dir is not None):
        for x in data[obsA].unique():
            fig, ax = plt.subplots(figsize=(10,10))
            
            df_plot = data[data[obsA]==x]
            
            if(obsB_colours is not None):
                obsB_colours = [obsB_colours[i] for i in df_plot[obsB]]
                
            ax.barh(df_plot[obsB], df_plot[obsB + "_frac"], color=obsB_colours)
            
            ax.set_title(obsA + ": " +  x + " (" + str(int(df_plot[obsA + "_ncells"].iloc[0]))+ " cells)")
            ax.set_xlabel("Fraction of cells")
            ax.set_ylabel(obsB)
            
            plt.tight_layout()
            
            # Replace / in obs names to avoid directory confusion
            x = x.replace('/', '_')
            
            plt.savefig(export_dir + x + ".pdf")
            plt.close(fig)
            
            

def makeAnnotationPlots(adata, clusters, model_predictions, markers, other_obs=None, export_dir= "", 
                        obsA_ignore_thresh = 50, obsB_ignore_thresh = 0.02):
    
    
    createAnnotationDirs(base_path = export_dir, clusters = clusters )
    
    
    # Plot UMAPs of SingleR predictions, clusterings
    sc.settings.figdir = export_dir
    sc.pl.umap(adata, color=[model_predictions] + other_obs + clusters,
               legend_loc="on data", legend_fontsize=4,
               save="_clusters.pdf",show=False)
    

    # Plot UMAPs of literature marker genes
    sc.settings.figdir = export_dir + "marker_expression/"
    for celltype in list(markers.keys()):
        sc.pl.umap(adata,color=markers[celltype],color_map="viridis_r",
                   save="_" + celltype.replace('/', '_') + "_markers.pdf",show=False)
    
    
    for cluster in clusters:
        
        # Plot bar charts of cluster prediction fractions
        plotMultiObsFraction(adata, obsA = cluster, obsB = model_predictions, 
                     obsA_ignore_thresh = 50 , obsB_ignore_thresh = 0.02,
                     #colours = scrabbit.getCelltypeColours(),
                     export_dir = export_dir + "cluster_fractions/" + cluster + "/")
        
        # Plot bar charts of prediction cluster fractions
        plotMultiObsFraction(adata, obsA = model_predictions, obsB = cluster, 
                     obsA_ignore_thresh = 50 , obsB_ignore_thresh = 0.02,
                     #colours = scrabbit.getCelltypeColours(),
                     export_dir = export_dir + "prediction_fractions/" + cluster + "/")
        


def filterMostCommonObs(adata,obs, frac_thresh=0.01):
    freq = adata.obs[obs].value_counts()
    freq_obs = freq[freq.gt(frac_thresh*freq.sum())].index.values
    return(freq_obs)




def assignCelltypeAdata(adata,clustering, clusters, slot_name, celltype_label):
    adata.obs.loc[adata.obs[clustering].isin(clusters),slot_name] = celltype_label
    

def assignCelltypeAdataSubset(adata, adata_sub, clustering, clusters, slot_name, celltype_label):
    adata.obs.loc[adata_sub.obs.index[adata_sub.obs[clustering].isin(clusters)],slot_name] = celltype_label
    assignCelltypeAdata(adata_sub, clustering, clusters, slot_name, celltype_label)
    
    



def plotObsUMAP3d(umap_3d, colour_by, hover_obs = None, layer=None,
                  colour_map=None,figsize=(1000,1000)):
    
    df = pd.DataFrame(umap_3d.obsm["X_umap"],columns=["UMAP1","UMAP2","UMAP3"])
    
    if((colour_by not in umap_3d.obs.columns) and (colour_by in umap_3d.var_names)):
        # TODO: Add continuous colour scale
        # TODO: Allow multiple genes
        df[colour_by] =  umap_3d.obs_vector(colour_by, layer=layer)
        exp_plot = True
    else:
        df[colour_by] = umap_3d.obs[colour_by].values
        exp_plot=False
        
    hover_str = None
    if(hover_obs is not None): 
        df[hover_obs] = umap_3d.obs[hover_obs].values
        hover_str = "".join(j + ": %{" + "customdata[" + 
                        "%d]} <br>" %k for k,j in enumerate(hover_obs)) + "<extra></extra>"        

    fig = go.Figure()
    
    if(exp_plot):
        fig = go.Figure(data=[go.Scatter3d(
            x=df["UMAP1"],
            y=df["UMAP2"],
            z=df["UMAP3"],
            mode='markers',
            customdata = df[hover_obs],
            marker=dict(
                size=2,
                color=df[colour_by],                # set color to an array/list of desired values
                colorscale='Viridis_r'),
            hovertemplate= hover_str
        )])
        
    else:
        colmap_dict=False
        colours = plot_utils.getClusterPalette()
        hover_data = None
        
        if(isinstance(colour_map,dict)):
           colmap_dict = True
           colours = colour_map
               
        for k,lbl in enumerate(df[colour_by].unique()):
            dfp = df[df[colour_by]==lbl]

            if(colmap_dict):
               col = colours[lbl]
            else:
               col = colours[k]
               
            if(hover_obs is not None):
                hover_data = dfp[hover_obs]

            fig.add_traces(go.Scatter3d(x=dfp['UMAP1'], y=dfp['UMAP2'],
                                        z=dfp["UMAP3"], mode='markers', 
                                        customdata = hover_data,
                                        name=lbl,
                                        marker = dict(color=col, size = 2),
                                        hovertemplate=hover_str))

    fig.update_layout(
        autosize=False,
        showlegend=True,
        legend= {'itemsizing': 'constant'},
        width=figsize[0],
        height=figsize[1],
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False)
        )
    )

    return(fig)
    