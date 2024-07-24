# Metabolights Metabolite-meta network

## Downloading studies
Study data is downloaded via FTP. We download the sample metadata file and the 'maf' files for all studies of interest. 

## Processing study data
Studies are objects of the `MTBLSDataset` class. This class has 3 functions:
- `read_data`
- `preprocess_data`
- `get_pathway_data`
- `plot_qc`
- `da_testing`

The data (maf files and metadata) are read in. If there is a single assay file, this is processed individually. If there are multiple maf files, they are processed invidiually, standardised, and concatenated at the end.

> If there are multiple maf files, the samples may have specific naming conventions e.g. 'X_POS' and 'X_NEG'. This needs to be harmonised into a unified identifier to continue with the integration. 

## Creating bipartite graph

## Visualising network
Network 