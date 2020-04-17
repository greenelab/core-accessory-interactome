# Using R 3.6.3
#install.packages("BiocManager")
#BiocManager::install("WGCNA")
#install.packages("flashClust")

# Load the WGCNA package
library(stats)
library(WGCNA)
library(flashClust)

get_threshold <- function(input_data_file){
    # Read in processed expression data

    # Note: We do not recommend attempting WGCNA on a data set consisting of fewer than 15 samples. 
    # In a typical high-throughput setting, correlations on fewer than 15 samples will simply be 
    # too noisy for the network to be biologically meaningful. If at all possible, one should have 
    # at least 20 samples; as with any analysis methods, more samples usually lead to more robust 
    # and refined results.

    # Read in data
    # Matrix that is samples x gene
    expression_data <- as.data.frame(read.table(input_data_file, sep="\t", header=TRUE, row.names=1))


    # Step 1: 
    # Start with similarity matrix, S = [s_ij] = |cor(i,j)| \in [0,1]
    # which represents the concordance between gene expression profiles for gene i and j.
    # S will eventually be transformed into adjacency matrix, A = [a_ij] 
    # which encodes the strength of the connection between gene i and j.
    # For 'hard thresholding' a_ij = 1 if s_ij >= threshold, else 0.
    # For 'soft thresholding' a_ij = |s_ij|^\beta .

    # We need to choose the soft thresholding power, \beta, to which
    # co-expression similarity is raised to calculate adjacency
    # Here, the soft thresholding power is chosen based on the criterion of 
    # approximate scale-free topology, pickSoftThreshold

    # Set of candidate powers
    powers = c(c(1:10), seq(from = 12, to=20, by=2))

    # Call the network topology analysis function
    # For each power the 'scale free topology fit index' is calculated 
    # and returned along with other information on connectivity.
    sft = pickSoftThreshold(expression_data, powerVector = powers, verbose = 5)

    # Plot the results:
    sizeGrWindow(9, 5)
    par(mfrow = c(1,2));
    cex1 = 0.9;

    # Plot soft-thresholding power parameter, \beta, wrt to scale-free topology fit index
    # Based on the plot, how to do we choose \beta where R^2 levels off and R^2 > 0.9

    # Criteria 1: Most biologists would be very suspicious of a gene co-expression network
    # that does not satisfy 'scale-free topology', at least approximately.
    # Scale-free topology is defined by a network where the probability that a node is 
    # connected with k other node (the degree distribution p(k) of a network) decays as a power law
    # p(k) ~ k^(-x).
    # Scale-free networks are extremely heterogeneous, their topology being dominated by a few highly
    # connected nodes (hubs) which link the rest of the less connected nodes to the system
    # Linear regression model fitting index R2 can be used to quantify how well a network 
    # satisfies a scale-free topology
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
        xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
        main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
        labels=powers,cex=cex1,col="red");
    abline(h=0.80,col="red")

    # Criteria 2: There is a natural trade-off between maximizing scale-free topology model
    # fit (R2) and maintaining a high mean number of connections: parameter values that 
    # lead to an R2 value close to 1 may lead to networks with very few connections
    # Plot mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
        xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
        main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

    # Save plots
}

generate_network_modules <- function(softPower,
input_data_file,
output_file){
    # Read in data
    # Matrix that is samples x gene
    expression_data <- as.data.frame(read.table(input_data_file, sep="\t", header=TRUE, row.names=1))

    # Step 2: 
    # Calculate adjacency matrix using selected soft-power
    #softPower = 20;
    adj= adjacency(expression_data,type = "unsigned", power = softPower)

    # Step 3: 
    # An important aim of co-expression network analysis is to detect subsets of
    # nodes (modules) that are tightly connected to each other
    # the adjacency matrix, A, is transformed into 
    # Topological overlap matrix (TOM) = [w_ij] = f(A) provides a similarity measure 
    # where dissimilarity = 1 - w_ij
    # This transformation to TOM is meant to minimize effects of noise and spurious associations, 
    # by taking into account neighboring genes (i.e. for each gene pair, looking at the similarity
    # profiles of their neighbors too)
    # The TOM combines the adjacency of two genes and the connection strengths these two genes share
    # with other "third party" genes
    # https://www.jstor.org/stable/pdf/2290430.pdf?casa_token=hWqqK4-65HYAAAAA:nNtlM8scsjjadud0O1eRw_TxIbgFGfr3tNr9Ev5COB2oEq3zcdEED0IhiOtbssfbtENeS1rT7TcYCZU7s2gythS7VGJ_oOQaoN7hL-b_le-lX3G216d_
    TOM = TOMsimilarity(adj);
    dissTOM = 1-TOM

    # Step 4: 
    # Use hierarchal clustering to produce a hierarchal clustering tree of genes (dendogram)
    # The dissimilarity measure is used to cluster genes using average linkage hierarchical clustering.
    # Hierarchal clustering merges closest clusters until there is only 1 cluster.
    # 'Average-link' clustering defines the distance between two clusters as the average
    # distance between their members.
    geneTree = hclust(as.dist(dissTOM), method = "average");
    # Plot the resulting clustering tree (dendrogram)
    sizeGrWindow(12,9)
    plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
        labels = FALSE, hang = 0.04);


    # Step 5:
    # Modules = branches of the hierarchal clustering tree
    minModuleSize = 30;
    # Module identification using dynamic tree cut:
    # WHAT IS THIS RETURNING????? <----------------------
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize);
    table(dynamicMods)

    # Plot modules on dendogram
    # Convert numeric lables into colors
    dynamicColors = labels2colors(dynamicMods)
    table(dynamicColors)
    # Plot the dendrogram and colors underneath
    sizeGrWindow(8,6)
    plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = "Gene dendrogram and module colors")


    # Step 6:
    # The Dynamic Tree Cut may identify modules whose expression profiles are very similar. 
    # It may be prudent to merge such modules since their genes are highly co-expressed. 
    # To quantify co-expression similarity of entire modules, we
    # calculate their eigengenes and cluster them on their correlation:

    # Calculate eigengenes
    MEList = moduleEigengenes(expression_data, colors = dynamicColors)
    MEs = MEList$eigengenes
    # Calculate dissimilarity of module eigengenes
    MEDiss = 1-cor(MEs);
    MEDissThres = 0.25
    # Plot the cut line into the dendrogram
    abline(h=MEDissThres, col = "red")
    # Call an automatic merging function
    merge = mergeCloseModules(expression_data, dynamicColors, cutHeight = MEDissThres, verbose = 3)
    # The merged module colors
    mergedColors = merge$colors;
    # Eigengenes of the new merged modules:
    mergedMEs = merge$newMEs;

    # Cluster module eigengenes
    METree = hclust(as.dist(MEDiss), method = "average");
    # Plot the result
    sizeGrWindow(7, 6)
    plot(METree, main = "Clustering of module eigengenes",
        xlab = "", sub = "")

    sizeGrWindow(12, 9)
    #pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
    plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                        c("Dynamic Tree Cut", "Merged dynamic"),
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)

    # Rename to moduleColors
    moduleColors = mergedColors
    # Construct numerical labels corresponding to the colors
    colorOrder = c("grey", standardColors(50));
    moduleLabels = match(moduleColors, colorOrder)-1;
    MEs = mergedMEs

    # Save module colors and labels for use in subsequent parts
    gene_ids <- colnames(expression_data)
    if (exists("mergedColors")){
    module_df <- as.data.frame(mergedColors, row.names=gene_ids)
    } else{
    module_df <- as.data.frame(dynamicColors, row.names=gene_ids)
    }
    write.table(module_df, output_file, append = FALSE, sep = "\t", 
                row.names = TRUE, col.names = TRUE)
}