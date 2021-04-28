# Network analysis
This directory generates G-G networks by applying clustering to the correlation matrix

[1_correlation_analysis.ipynb](1_correlation_analysis.ipynb) Performs correlation analysis to compare the similarity between genes and applies different threshold cutoffs to determine the strength of connection between genes. This notebook also visualizes the correlation structure to inform the clustering approach used in the next notebook.

[2_get_network_communities.ipynb](2_get_network_communities.ipynb) Applies clustering methods to identify gene modules based on the correlation matrices.

[3_module_validation.ipynb](3_module_validation.ipynb) Validates generated modules.

<!---Say which clustering method was decided and why>
<!---Plots of validation>