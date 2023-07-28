# 2013TSAGWAS-PheWASResultsViewer.py
# PheWAS and IUT Results Viewer of 2013TSA GWAS Hits
# Alexander Eyre (2023)

#### Define Functions ####

def manhattan(df, pvalues, traits, colors, curr_chr, bct05, bct10, dataset):

    fig           = matplotlib.pyplot.figure(figsize=(18, 4))
    
    x_vals        = list(np.arange(1,len(df)+1))
    y_vals        = -np.log10(list(df[pvalues]))
    
    x_plot        = np.linspace(0, 100, 1000)
    y_plot05      = [i*-np.log10(bct05) for i in np.ones(1000)]
    #y_plot10      = [i*-np.log10(bct10) for i in np.ones(1000)]

    ax            = sns.scatterplot(x = x_vals, 
                                    y = y_vals, 
                                    hue = df[colors], 
                                    palette = 'tab20', 
                                    edgecolor = "#000000",
                                    linewidth = 0.5)
    
    ax.plot(x_plot, y_plot05, color = 'r', linestyle = 'dashed')
    #ax.plot(x_plot, y_plot10, color = 'b', linestyle = 'dashed')
    
    ax.set_title("%s" % dataset, fontsize = 14)
    
    ax.set_ylabel("-log10(p)", fontsize = 12)
    ax.set_yticks([0,1,2,3,4,5])
    
    ax.set_xticks(x_vals)
    ax.tick_params(bottom=True)
    ax.tick_params(axis='x', colors='#D0D0D0')
    ax.set_xticklabels(list(df[traits]), fontsize = 8, rotation = 60, ha = 'right', color = 'black')
    ax.set_xlabel("CBARQ Questions", fontsize = 15)
    ax.set_xlim(left = 0, right = len(df)+1)
    
    handles, labels = ax.get_legend_handles_labels()
    for ha in handles:
        ha.set_edgecolor("#000000")
        ha.set_linewidth(0.5)
    
    ax.legend(handles, labels, loc = (1.0,-0.2), labelspacing = 1.0, frameon = False, edgecolor = "black").remove()
    
def volcano(df, traits, pvalues, betavalues, colors, curr_chr, bct05, dataset):
    
    fig           = matplotlib.pyplot.figure(figsize=(8, 6))

    ax            = sns.scatterplot(x = list(df[betavalues]), 
                                    y = list(df[pvalues]), 
                                    hue = df[colors], 
                                    palette = 'tab20', 
                                    edgecolor = "#000000", 
                                    linewidth = 0.5,
                                    style = df[colors])
    
    stepsize      = (max(df[betavalues]) - min(df[betavalues])) * 0.1
    
    x_plot        = np.linspace(min(df[betavalues]) - stepsize, max(df[betavalues]) + stepsize, 1000)
    y_plot05      = [i*-np.log10(bct05) for i in np.ones(1000)]
    y_plot0       = np.arange(-0.2, 7.2, (7.2+0.2)/999.5)
    
    i = 0
    while i < len(list(df[betavalues])):
        if list(df[pvalues])[i] > -np.log10(bct05):
            ax.text(list(df[betavalues])[i] + stepsize * 0.1, 
                    list(df[pvalues])[i],
                    list(df[traits])[i], 
                    horizontalalignment = 'left',
                    size = 'small',
                    color = 'black')
        i = i + 1
    
    ax.plot(np.zeros(1000), y_plot0, color = 'gray')
    ax.plot(x_plot, y_plot05, color = 'r', linestyle = 'dashed')
    
    ax.set_title("%s" % dataset, fontsize = 14)
    
    ax.set_ylabel("-log10(p)", fontsize = 12)
    ax.set_yticks([0,1,2,3,4,5,6,7])
    ax.set_ylim(-0.1, np.ceil(max(df[pvalues])) + 0.1)
    
    ax.set_xlabel("Beta", fontsize = 12)
    ax.set_xlim(left = min(df[betavalues])-stepsize, right = max(df[betavalues])+stepsize)

    handles, labels = ax.get_legend_handles_labels()
    for ha in handles:
        ha.set_edgecolor("#000000")
        ha.set_linewidth(0.5)
    
    ax.legend(handles, labels, loc = (1.0,0.1), labelspacing = 1.0, frameon = False, edgecolor = "black").remove()
    
def compareplots(genotype_data, CBARQ_data, chr_name, CBARQ_Q, source):    
    
    fig            = matplotlib.pyplot.figure(figsize=(8, 6))
    
    # Identify all of the unique genotypes
    genotype_list  = []
    for x in genotype_data[chr_name]:
        if x not in genotype_list:
            genotype_list.append(x)
    
    genotype_list.sort()
    genotype_count = len(genotype_list)    
    
    # Generate Datasets
    i              = 0
    dataset        = {}
    curr_vals      = []
    while i < genotype_count:
        curr_vals = list(CBARQ_data[CBARQ_Q][genotype_data[chr_name] == genotype_list[i]])
        if "nan" in curr_vals:
            curr_vals.remove("nan")
        curr_vals = [x for x in curr_vals if np.isnan(x) == False]
        curr_vals = [x + i / genotype_count for x in curr_vals]   
        dataset[genotype_list[i]] = curr_vals
        i = i + 1  
    
    # Remove empty datasets
    for genotype in genotype_list:
        if dataset[genotype] == []:
            del dataset[genotype]
            genotype_list.remove(genotype)
            genotype_count = genotype_count - 1
    
    # Generate Plot
    i              = 0
    while i < genotype_count:
        hplot = sns.histplot(data    = dataset[genotype_list[i]],
                             bins    = np.arange(1, 6.5, 1 / genotype_count),
                             alpha   = 0.8,
                             stat    = "probability") 
        i = i + 1
    
    # Generate Legend
    legend_artist  = [Rectangle((0, 0), 1, 1, fc = '#1f77b4', alpha = 0.8),
                      Rectangle((0, 0), 1, 1, fc = '#ff7f0e', alpha = 0.8),
                      Rectangle((0, 0), 1, 1, fc = '#2ca02c', alpha = 0.8),
                      Rectangle((0, 0), 1, 1, fc = '#d62728', alpha = 0.8)] 
                      
    if genotype_count == 2:
        matplotlib.pyplot.legend(handles = [legend_artist[0], legend_artist[1]],
                                 labels  = ["%s (%s)" % (genotype_list[0], len(dataset[genotype_list[0]])), 
                                            "%s (%s)" % (genotype_list[1], len(dataset[genotype_list[1]]))])
    elif genotype_count == 3:
        matplotlib.pyplot.legend(handles = [legend_artist[0], legend_artist[1], legend_artist[2]],
                                 labels  = ["%s (%s)" % (genotype_list[0], len(dataset[genotype_list[0]])), 
                                            "%s (%s)" % (genotype_list[1], len(dataset[genotype_list[1]])),
                                            "%s (%s)" % (genotype_list[2], len(dataset[genotype_list[2]]))])
    else:
        matplotlib.pyplot.legend(handles = [legend_artist[0], legend_artist[1], legend_artist[2], legend_artist[3]],
                                 labels  = ["%s (%s)" % (genotype_list[0], len(dataset[genotype_list[0]])), 
                                            "%s (%s)" % (genotype_list[1], len(dataset[genotype_list[1]])),
                                            "%s (%s)" % (genotype_list[2], len(dataset[genotype_list[2]])),
                                            "%s (%s)" % (genotype_list[3], len(dataset[genotype_list[3]]))])
    # Additional Plot Settings
    matplotlib.pyplot.title(source, size = 14)
    matplotlib.pyplot.xticks(ticks  = [1.5, 2.5, 3.5, 4.5, 5.5],
               labels = ["1", "2", "3", "4", "5"])
    matplotlib.pyplot.xlabel("Score", fontsize = 12)
    matplotlib.pyplot.ylabel("Relative Frequency", fontsize = 12)

def compareplotsAUS(genotype_data, CBARQ_data, chr_name, CBARQ_Q, source):    
    
    fig            = matplotlib.pyplot.figure(figsize=(8, 6))
    
    # Identify all of the unique genotypes
    genotype_list  = []
    for x in genotype_data[chr_name]:
        if x not in genotype_list:
            genotype_list.append(x)
    
    genotype_list.sort()
    genotype_count = len(genotype_list)    
    
    # Generate Datasets
    i              = 0
    dataset        = {}
    curr_vals      = []
    while i < genotype_count:
        curr_vals = list(CBARQ_data[CBARQ_Q][genotype_data[chr_name] == genotype_list[i]])
        if "nan" in curr_vals:
            curr_vals.remove("nan")
        curr_vals = [x for x in curr_vals if np.isnan(x) == False]
        curr_vals = [x + i / genotype_count for x in curr_vals]   
        dataset[genotype_list[i]] = curr_vals
        i = i + 1  
    
    # Remove empty datasets
    for genotype in genotype_list:
        if dataset[genotype] == []:
            del dataset[genotype]
            genotype_list.remove(genotype)
            genotype_count = genotype_count - 1
    
    # Generate Plot
    i              = 0
    while i < genotype_count:
        hplot = sns.histplot(data    = dataset[genotype_list[i]],
                             bins    = np.arange(0, 5.5, 1 / genotype_count),
                             alpha   = 0.8,
                             stat    = "probability") 
        i = i + 1
    
    # Generate Legend
    legend_artist  = [Rectangle((0, 0), 1, 1, fc = '#1f77b4', alpha = 0.8),
                      Rectangle((0, 0), 1, 1, fc = '#ff7f0e', alpha = 0.8),
                      Rectangle((0, 0), 1, 1, fc = '#2ca02c', alpha = 0.8),
                      Rectangle((0, 0), 1, 1, fc = '#d62728', alpha = 0.8)] 
                      
    if genotype_count == 2:
        matplotlib.pyplot.legend(handles = [legend_artist[0], legend_artist[1]],
                                 labels  = ["%s (%s)" % (genotype_list[0], len(dataset[genotype_list[0]])), 
                                            "%s (%s)" % (genotype_list[1], len(dataset[genotype_list[1]]))])
    elif genotype_count == 3:
        matplotlib.pyplot.legend(handles = [legend_artist[0], legend_artist[1], legend_artist[2]],
                                 labels  = ["%s (%s)" % (genotype_list[0], len(dataset[genotype_list[0]])), 
                                            "%s (%s)" % (genotype_list[1], len(dataset[genotype_list[1]])),
                                            "%s (%s)" % (genotype_list[2], len(dataset[genotype_list[2]]))])
    else:
        matplotlib.pyplot.legend(handles = [legend_artist[0], legend_artist[1], legend_artist[2], legend_artist[3]],
                                 labels  = ["%s (%s)" % (genotype_list[0], len(dataset[genotype_list[0]])), 
                                            "%s (%s)" % (genotype_list[1], len(dataset[genotype_list[1]])),
                                            "%s (%s)" % (genotype_list[2], len(dataset[genotype_list[2]])),
                                            "%s (%s)" % (genotype_list[3], len(dataset[genotype_list[3]]))])
    # Additional Plot Settings
    matplotlib.pyplot.title(source, size = 14)
    matplotlib.pyplot.xticks(ticks  = [0.5, 1.5, 2.5, 3.5, 4.5],
               labels = ["1", "2", "3", "4", "5"])
    matplotlib.pyplot.xlabel("Score", fontsize = 12)
    matplotlib.pyplot.ylabel("Relative Frequency", fontsize = 12)

#### Import Packages ####

import pandas as pd
import numpy as np
import matplotlib
import seaborn as sns
import streamlit as st
from matplotlib.patches import Rectangle
sns.set(style = "whitegrid")

#### Data Import & Processing ####

# Import PheWAS Results for each dataset
AUS_data          = pd.read_excel("PheWAS Results.xlsx", sheet_name = "AUS")
UK_data           = pd.read_excel("PheWAS Results.xlsx", sheet_name = "UK")
SE_data           = pd.read_excel("PheWAS Results.xlsx", sheet_name = "SE")

# Import IUT Results for each Fisher's combined dataset
AUSUK_data        = pd.read_excel("PheWAS Results.xlsx", sheet_name = "AUS-UK")
AUSSE_data        = pd.read_excel("PheWAS Results.xlsx", sheet_name = "AUS-SE")
UKSE_data         = pd.read_excel("PheWAS Results.xlsx", sheet_name = "UK-SE")
AUSUKSE_data      = pd.read_excel("PheWAS Results.xlsx", sheet_name = "AUS-UK-SE")

# Import Bonferroni Corrected Thresholds for each Fisher's combined dataset
AUS_BCT           = list(AUS_data["Bonferroni Corrected Threshold"].dropna())
UK_BCT            = list(UK_data["Bonferroni Corrected Threshold"].dropna())
SE_BCT            = list(SE_data["Bonferroni Corrected Threshold"].dropna())
AUSUK_BCT         = list(AUSUK_data["Bonferroni Corrected Threshold"].dropna())
AUSSE_BCT         = list(AUSSE_data["Bonferroni Corrected Threshold"].dropna())
UKSE_BCT          = list(UKSE_data["Bonferroni Corrected Threshold"].dropna())
AUSUKSE_BCT       = list(AUSUKSE_data["Bonferroni Corrected Threshold"].dropna())

# Import CBARQ Factors 
CBARQ_factors     = pd.read_excel("PheWAS Results.xlsx", sheet_name = "CBARQ Qs")
CBARQ_factors     = CBARQ_factors.iloc[:,[0,2]]
CBARQ_factors.index = list(np.arange(1,101,1))

# Generate an ordering system for CBARQ Factors
trait_order       = {"Trainability":1, "OwnDirAgg":2, "StrDirAgg":3, "DogDirAgg":4,
                     "FamDogAgg":5, "StrDirFear":6, "NonSocFear":7, "DogDirFear":8,
                     "SepRelBeh":9, "AtcAtnSeek":10, "TouchSens":11, "Chasing":12,
                     "Excitability":13, "Misc":14}
trait_list        = ["Trainability", "OwnDirAgg", "StrDirAgg", "DogDirAgg", "FamDogAgg",  
                     "StrDirFear", "NonSocFear", "DogDirFear", "SepRelBeh", "AtcAtnSeek",
                     "TouchSens", "Chasing", "Excitability",  "Misc"]
trait_order       = {nm:i for i,nm in enumerate(trait_list)}

# Import CTC Factors
CTC_factors       = pd.read_excel("TSA Pretraining Factors.xlsx", sheet_name = "Sheet1")
CTC_questions     = CTC_factors.iloc[:,5]
CTC_factors       = CTC_factors.iloc[:,[0,2]]
CTC_factors.index = list(np.arange(1,101,1))

# Generate an ordering system for CTC Factors
CTC_order         = {"Trainability":0, "Handler Aggression":1, "Stranger Aggression":2,
                    "Dog Aggression":3, "Stranger Fear":4, "Environmental Fear":5,
                    "Dog Fear": 6, "Touch Sensitivity":7, "Separation Behavior":8, 
                    "Excitability":9, "Attention Seeking":10, "Distraction":11, "Misc":12}
CTC_list          = ["Trainability", "Handler Aggression", "Stranger Aggression",
                    "Dog Aggression", "Stranger Fear", "Environmental Fear",
                    "Dog Fear", "Touch Sensitivity", "Separation Behavior", 
                    "Excitability", "Attention Seeking", "Distraction", "Misc"]
CTC_order         = {nm:i for i, nm in enumerate(CTC_list)}

# Import Alleles & Genotypes
AUS_alleles       = pd.read_table("AUSAlleles.txt")
AUS_genotypes     = pd.read_table("AUSGenotypes.txt")
SE_alleles        = pd.read_table("SEAlleles.txt")
SE_genotypes      = pd.read_table("SEGenotypes.txt")
UK_alleles        = pd.read_table("UKAlleles.txt")
UK_genotypes      = pd.read_table("UKGenotypes.txt")

# Import CBARQ Scores
AUS_CBARQ         = pd.read_table("AUS_PostCBARQ.txt", sep = ",")
UK_CBARQ          = pd.read_table("UK_PostCBARQ.txt", sep = ",")
SE_CBARQ          = pd.read_table("SE_PostCBARQ.txt", sep = ",")

# Generate SNP and CBARQ Question Lists
SNP_list          = ("chr1.22989459", "chr1.24927539", "chr1.25289424", "chr3.47103534", "chr3.47134935", "chr3.47212502",
                     "chr4.31420247", "chr6.76632282", "chr7.66358701", "chr9.16559175", "chr11.38609196",
                     "chr13.55534649", "chr13.58498102", "chr13.59658137", "chr13.59763520", "chr13.59902870", 
                     "chr15.40757218", "chr19.21040815", "chr36.25252101", "chr36.25648690")
CBARQ_list        = list(AUS_CBARQ.columns)[1:101]

#### Steamlit App ####
## Compile the sidebar for SNP and CBARQ selection
st.set_option('deprecation.showPyplotGlobalUse', False)
st.set_page_config(layout="wide")
base = "dark"
st.subheader("PheWAS & IUT Results Viewer of 2013TSA GWAS Hits")
st.caption("Alexander Eyre (2023)")
st.sidebar.markdown("Select SNP and CBARQ Question")

# Choose SNP to view
SNP_selector      = st.sidebar.selectbox("SNP", SNP_list)

# Choose CBARQ Question to view
CBARQ_selector    = st.sidebar.selectbox("CBARQ Question", CBARQ_list)

# Explanation for the CBARQ Question
st.sidebar.markdown("%s" % list(CTC_factors["TSA_Pretraining_Factor"])[list(CTC_factors["Question#"]).index(CBARQ_selector)])
st.sidebar.markdown("%s" % list(CTC_questions)[list(CTC_factors["Question#"]).index(CBARQ_selector)])

# Add Legend for Figures
st.sidebar.image("FigLegend.jpg", caption = "Legend")

## Generate the webapp
volc, histo        = st.columns(2)

if SNP_selector in AUS_data.columns:
    if CBARQ_selector in AUS_CBARQ.columns:
        AUS_snp_pos    = list(AUS_data.columns).index(SNP_selector)
        AUS_df         = AUS_data.iloc[1:,[0, AUS_snp_pos+1, AUS_snp_pos+2]]
        AUS_df.columns = ["Phenotype", "-Log(P)", "Beta"]
        AUS_df         = AUS_df.join(CTC_factors["TSA_Pretraining_Factor"])
        AUS_df         = AUS_df.sort_values(by = 'TSA_Pretraining_Factor', key = lambda nm: nm.map(CTC_order))
        AUS_df         = AUS_df.dropna()
        volc.pyplot(fig = volcano(AUS_df, "Phenotype", "-Log(P)", "Beta", "TSA_Pretraining_Factor", SNP_selector, AUS_BCT[0], "Australian"),
                  clear_figure = True)
              
if SNP_selector in AUS_data.columns:
    if CBARQ_selector in AUS_CBARQ.columns:
        histo.pyplot(fig = compareplotsAUS(AUS_genotypes, AUS_CBARQ, SNP_selector, CBARQ_selector, "Australian"),
                  clear_figure = True)
                   
if SNP_selector in UK_data.columns:
    if CBARQ_selector in UK_CBARQ.columns:
        UK_snp_pos     = list(UK_data.columns).index(SNP_selector)
        UK_df          = UK_data.iloc[1:,[0, UK_snp_pos+1, UK_snp_pos+2]]
        UK_df.columns  = ["Phenotype", "-Log(P)", "Beta"]
        UK_df          = UK_df.join(CTC_factors["TSA_Pretraining_Factor"])
        UK_df          = UK_df.sort_values(by = 'TSA_Pretraining_Factor', key = lambda nm: nm.map(CTC_order))
        UK_df          = UK_df.dropna()
        volc.pyplot(fig = volcano(UK_df, "Phenotype", "-Log(P)", "Beta", "TSA_Pretraining_Factor", SNP_selector, UK_BCT[0], "United Kingdom"),
                  clear_figure = True)
              
if SNP_selector in UK_data.columns:
    if CBARQ_selector in UK_CBARQ.columns:
        histo.pyplot(fig = compareplots(UK_genotypes, UK_CBARQ, SNP_selector, CBARQ_selector, "United Kingdom"),
                  clear_figure = True)  
              

if SNP_selector in SE_data.columns:
    if CBARQ_selector in SE_CBARQ.columns:
        SE_snp_pos     = list(SE_data.columns).index(SNP_selector)
        SE_df          = SE_data.iloc[1:,[0, SE_snp_pos+1, SE_snp_pos+2]]
        SE_df.columns  = ["Phenotype", "-Log(P)", "Beta"]
        SE_df          = SE_df.join(CTC_factors["TSA_Pretraining_Factor"])
        SE_df          = SE_df.sort_values(by = 'TSA_Pretraining_Factor', key = lambda nm: nm.map(CTC_order))
        SE_df          = SE_df.dropna()
        volc.pyplot(fig = volcano(SE_df, "Phenotype", "-Log(P)", "Beta", "TSA_Pretraining_Factor", SNP_selector, SE_BCT[0], "Seeing Eye"),
                  clear_figure = True)
                  
if SNP_selector in SE_data.columns:
    if CBARQ_selector in SE_CBARQ.columns:
        histo.pyplot(fig = compareplots(SE_genotypes, SE_CBARQ, SNP_selector, CBARQ_selector, "Seeing Eye"),
                  clear_figure = True)
                  
# AUS & UK
if SNP_selector in AUS_data.columns and SNP_selector in UK_data.columns:
    AUSUK_df       = AUSUK_data.iloc[1:,[0, list(AUSUK_data.columns).index(SNP_selector) + 3]]
    AUSUK_df.columns = ["Phenotype", "-Log(P)"]
    AUSUK_df       = AUSUK_df.join(CTC_factors["TSA_Pretraining_Factor"])
    AUSUK_df       = AUSUK_df.dropna()  
    AUSUK_df       = AUSUK_df.sort_values(by = 'TSA_Pretraining_Factor', key = lambda nm: nm.map(CTC_order)) 
    st.pyplot(fig = manhattan(AUSUK_df, "-Log(P)", "Phenotype", "TSA_Pretraining_Factor", SNP_selector, AUSUK_BCT[0], AUSUK_BCT[1], "Australian - United Kingdom"),
                  clear_figure = True)
    
# AUS & SE
if SNP_selector in AUS_data.columns and SNP_selector in SE_data.columns:
    AUSSE_df       = AUSSE_data.iloc[1:,[0, list(AUSSE_data.columns).index(SNP_selector) + 3]]
    AUSSE_df.columns = ["Phenotype", "-Log(P)"]
    AUSSE_df       = AUSSE_df.join(CTC_factors["TSA_Pretraining_Factor"])
    AUSSE_df       = AUSSE_df.dropna() 
    AUSSE_df       = AUSSE_df.sort_values(by = 'TSA_Pretraining_Factor', key = lambda nm: nm.map(CTC_order)) 
    st.pyplot(fig = manhattan(AUSSE_df, "-Log(P)", "Phenotype", "TSA_Pretraining_Factor", SNP_selector, AUSSE_BCT[0], AUSSE_BCT[1], "Australian - Seeing Eye"),
                  clear_figure = True)    
# UK & SE
if SNP_selector in UK_data.columns and SNP_selector in SE_data.columns:
    UKSE_df        = UKSE_data.iloc[1:,[0, list(UKSE_data.columns).index(SNP_selector) + 3]]
    UKSE_df.columns = ["Phenotype", "-Log(P)"]
    UKSE_df        = UKSE_df.join(CTC_factors["TSA_Pretraining_Factor"])
    UKSE_df        = UKSE_df.dropna() 
    UKSE_df        = UKSE_df.sort_values(by = 'TSA_Pretraining_Factor', key = lambda nm: nm.map(CTC_order)) 
    st.pyplot(fig = manhattan(UKSE_df, "-Log(P)", "Phenotype", "TSA_Pretraining_Factor", SNP_selector, UKSE_BCT[0], UKSE_BCT[1], "United Kingdom - Seeing Eye"),
                  clear_figure = True) 
# AUS & UK & SE    
if SNP_selector in AUS_data.columns and SNP_selector in UK_data.columns and SNP_selector in SE_data.columns:   
    AUSUKSE_df     = AUSUKSE_data.iloc[1:,[0, list(AUSUKSE_data.columns).index(SNP_selector) + 5]]
    AUSUKSE_df.columns = ["Phenotype", "-Log(P)"]
    AUSUKSE_df     = AUSUKSE_df.join(CTC_factors["TSA_Pretraining_Factor"])
    AUSUKSE_df     = AUSUKSE_df.dropna()
    AUSUKSE_df     = AUSUKSE_df.sort_values(by = 'TSA_Pretraining_Factor', key = lambda nm: nm.map(CTC_order)) 
    st.pyplot(fig = manhattan(AUSUKSE_df, "-Log(P)", "Phenotype", "TSA_Pretraining_Factor", SNP_selector, AUSUKSE_BCT[0], AUSUKSE_BCT[1], "Australian - United Kingdom - Seeing Eye"),
                  clear_figure = True) 
