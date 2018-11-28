#!/bin/bash

PATH=$PATH 

export PATH 

computeMatrix scale-regions -S *.bw -R lhasa2_specific_cmt2-3.fx.DMR.bed lhasa2_olp_cmt2-3.fx.DMR.bed cmt2-3_specific_lhasa2.fx.DMR.bed --binSize 10 -b 1000 -a 1000 -p 3 --verbose --regionBodyLength 200 -out matrix_200.mat.gz

# -R   # regions to plot over -R
#-regionbodylength   # controls the length of the region you want all of the bed regions to be scaled to.
#-bin  # drastically changes the amount of noise, if the plots look very noisy, increase bin size
#-S   # list of single direction (+ strand converted DNA methylation) .BigWig (.bw) files in the folder you want plotted over the regions 
#-b  # area b4 the gene/region to derive counts from
#-a  # area after the gene/region to derive counts from

# Profile of above data
plotProfile -m matrix_200.mat.gz --legendLocation upper-right --outFileName Lhasa2_DMR_plot_200.eps --perGroup  --plotFileFormat eps --dpi 400  

#--perGroup   # splits the plots for each bed region into their own graphs. To plot all on the same (not always advisable biologically, bc your not comparing the same places) 
#--plotTitle   # "Hypomethylated_DMRs"
#--colors    # red yellow blue
#--plotType=fill  # add color between the x axis and the lines

# heatmap of same data
plotProfile -m matrix_200.mat.gz --perGroup --kmeans 2 --plotType heatmap --plotFileFormat eps --dpi 400 --outFileName Lhasa2_DMR_heatmap_200.eps 

