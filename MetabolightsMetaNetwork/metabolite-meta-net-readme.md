# Metabolights metabolite-meta network
Last updated 24/07/24

COVID-19 example notebook is stored in `MetHDN/MetDMN/COVID_network-cleaned.ipynb`

## Downloading studies
Study data is downloaded via FTP. We download the sample metadata file and the 'maf' files for all studies of interest. 

## Processing study data
Studies are objects of the `MTBLSDataset` class. This class has several functions:
- `read_data`
- `preprocess_data`
- `get_pathway_data`
- `plot_qc`
- `da_testing`

The data (maf files and metadata) are read in. If there is a single assay file, this is processed individually. If there are multiple maf files, they are processed invidiually, standardised, and concatenated at the end.

### Class attributes 
For each study the following is specified:
- File path containing maf files and metadata
- Metabolights ID
- Metadata group - column name to use in metadata file
- Metadata filter - keep only samples of the specified metadata labels (e.g. only 'Case' and 'Control' but not 'QC'). The user has to specify which label corresponds to the 'Case' and which the 'Control' so these can be standardised across studies. 
- Outliers (optional): outlier sample IDs to filter out

> If there are multiple maf files, the samples may have specific naming conventions e.g. 'X_POS' and 'X_NEG'. This needs to be harmonised into a unified identifier to continue with the integration. For now, the `remove_suffix` attribute has been added to the `MTBLSDataset` class where the user can specify how many characters need removing from the end of the sample name strings in both maf files and metadata. 

### Preprocessing and QC steps
Takes place within the `preprocess_data` function:

1. ChEBI `database_identifier` column is set as index. Any rows without a ChEBI are dropped.
2. Keep only columns (samples) present in the metadata filter 
3. Transpose the DataFrame such that rows represent samples and columns represent metabolites
4. Drop any outlier samples
5. Replace empty strings with `np.nan`
6. Drop rows/columns where all values are 0
7. Drop columns where over 50% of the samples have missing data
8. Compute/output level of missingness (%)
9. Impute missing values using `scikit-learn KNNImputer` 
10. Log2 transform imputed data
11. Auto-scale data with `scikit-learn StandardScaler`
12. Add `Group` metadata column at the end

Return processed DataFrame (`self.processed_data`)

QC can be checked at this point - the `plot_qc()` function outputs a basic PCA biplot coloured by the metadata grouping, as well as a boxplot of every n^th metabolite. 

> Further quality control steps may be recommended such as checking the PCA residuals, etc.

### Storing a pathway-space representation of the data
Reactome pathways are used to create single-sample pathway scores (ssPA) for the processed metabolite data. The `sspa` Python package we developed is used for this. Some pathway coverage statistics are also stored. 

A matrix of samples-by-pathways is returned and stored as a class attribute (`self.pathway_data`). 

### Differential abundance testing (metabolite or pathway level)
Either metabolite abundances or pathway scores can be tested for differential abundance (case vs. control groups). 

## Creating bipartite graph
The bipartite graph is created using the `get_bipartite()` function. This function takes a single argument `studies`: a list of studies. All objects of the `MTBLSDataset` class have several attributes including `node_name` which is the MetaboLights identifier, `DA_metabolites` which is a list of all the differential metabolites for the study and `connection` which represents a list of tuples which show the connection between each study and its differential metabolites. These are used to build a networkX bipartite graph. 

A bi-adjacency matrix can be created from the bipartite graph: this is a binary matrix represented as a DataFrame where columns represent differential metabolites and rows represent studies. The presence of a 1 or 0 indicates whether a metabolite is differential in a study. 

The adjacency matrix can be created by multiplying the by transposed bi-adjacency matrix by itself (biadj.T*biadj). This matrix shows the co-occurrence of differential metabolites across studies. An edge list DataFrame is created from the adjacency matrix: it contains 3 columns: `source`, `target`, and `weight`. The `weight` represents the number of studies a metabolite co-occurs in. Next the metabolite meta-network is created using the networkX `from_pandas_edgelist()` function. Self loops are removed. Node attributes are added, for example ChEBI name, study contributions, and direction of metabolite differential abundance. This graph can then be exported as a `.graphML` file or other formats for further visualisation for example Cytoscape.


## Visualising network
NetworkX is exported to Cytoscape - ensure key attributes such as study contribution, ChEBI name, etc are added to network (see COVID notebook example).