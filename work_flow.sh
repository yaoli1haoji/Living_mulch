#!/bin/bash

#import the raw data. 
#They are already demutiplexed Pairedend data
#Use According import format


#import the raw data. 
#They are already demutiplexed Pairedend data
#Use According import format

#set_up_raw_fasta_file_directory 
fastqDir=living_mulch_data

#set_up_trimmed__fasta_file_directory
trimDir=trimmed_data
absolye_path=/home/hl46161/publish_living_mulch/trimmed_data


#First thing to do is to trim the primer sequence in the seuquence 
for sample in $(ls $fastqDir| grep .*R1_001.fastq.gz | awk -F_ '{print $1"_"$2}'); do
       
       echo $sample is being processed
       
       cutadapt -q 20 -m 1 -g CTGTCTCTTATACACATCT -G GAGGACTACNVGGGTWTCTAAT -o $trimDir/${sample}_L001_R1_001.fastq.gz -p $trimDir/${sample}_L001_R2_001.fastq.gz $fastqDir/${sample}_L001_R1_001.fastq.gz $fastqDir/${sample}_L001_R2_001.fastq.gz

done

### create manifest_file that are required for qiime2 data import
 for sample in $(ls $trimDir); do
       
       echo $sample
       if [[ $sample == *_"R1"_* ]];
       then
          revised_sample=`echo $sample | cut -d'_' -f 1`
          echo $revised_sample
          echo $revised_sample,$absolye_path/$sample,forward >> manifest_file_test.txt
       else
          revised_sample=`echo $sample | cut -d'_' -f 1`
          echo $revised_sample
          echo $revised_sample,$absolye_path/$sample,reverse >> manifest_file_test.txt
       fi
done


#import the 16s illuma reads data as paired end format 

 qiime tools import \
     --type 'SampleData[PairedEndSequencesWithQuality]' \
     --input-path manifest_file.csv \
     --output-path paired-end-demux.qza \
     --input-format PairedEndFastqManifestPhred33

#produce a visualization file 

 qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization demux.qzv

#
 #join the paired reads together 
    qiime vsearch join-pairs --i-demultiplexed-seqs paired-end-demux.qza --o-joined-sequences living_mulch_trimmed_jointed_data.qza
    
         qiime demux summarize \
      --i-data living_mulch_trimmed_jointed_data.qza \
      --o-visualization living_mulch_trimmed_jointed_data.qzv


######################################################
      
# # # # #prefilter the reads with quality score
   qiime quality-filter q-score \
       --i-demux living_mulch_trimmed_jointed_data.qza \
       --o-filtered-sequences living_mulch_trimmed_jointed_data_filtered_data.qza \
       --o-filter-stats living_mulch_trimmed_jointed_data_filtered_data_stats.qza
       
# # vivualize the result after quality filtering      
   qiime metadata tabulate \
    --m-input-file living_mulch_trimmed_jointed_data_filtered_data_stats.qza \
    --o-visualization living_mulch_trimmed_jointed_data_filtered_data_stats.qzv
      
# # #generate the feature table 
   qiime deblur denoise-16S \
      --i-demultiplexed-seqs living_mulch_trimmed_jointed_data_filtered_data.qza \
      --p-trim-length 240 \
      --p-sample-stats \
      --o-representative-sequences rep-seqs.qza \
      --o-table table.qza \
      --o-stats deblur-stats.qza

      # # summarize the original feature table
    qiime feature-table summarize \
       --i-table table.qza \
       --o-visualization table.qzv\
    
       # summarize the original represntative sequence   
    qiime feature-table tabulate-seqs \
      --i-data rep-seqs.qza \
      --o-visualization rep-seqs.qzv
      
    #generate the phylogenetic tree 
     qiime phylogeny align-to-tree-mafft-fasttree \
     --i-sequences rep-seqs.qza \
     --o-alignment aligned-rep-seqs.qza \
     --o-masked-alignment masked-aligned-rep-seqs.qza \
     --o-tree unrooted-tree.qza \
     --o-rooted-tree rooted-tree.qza
     

     #####alpha and beta diversity test of original table at 20000 seuqence 
        qiime diversity core-metrics-phylogenetic \
   --i-phylogeny rooted-tree.qza \
   --i-table table.qza \
   --p-sampling-depth 20000 \
   --m-metadata-file living_mulch_data_metafile_massalin.tsv \
   --output-dir core-metrics-results-orginal-table-20000
   
   ## according to the visualization of beta diveristy of samples, three blank samples contain higher sequences 
   ## therefore I need to decontam the data 
   ## first thing I did is to remove all the sequences of P2BI from the whole dataset
   ## since P2BI located far away from all the other samples based on beta diveristy, filtered out its sequences will not make big change to the datasets. (I already test it)
   
  #exported table.qza for filtering contanmination
  qiime tools export \
    --input-path table.qza \
    --output-path exported-feature-table
    
  #convert biom file to txt file, a format which readbale for script in pycharm 
  biom convert -i feature-table.biom -o original_table.txt --to-tsv
     

  #after the filtration is done in pycharm, convert the filtered feature table from txt to biom blank sample and inrow sample was filtered 
   biom convert -i /home/hl46161/publish_living_mulch/exported-feature-table/table_P2BI_filtered.txt -o /home/hl46161/publish_living_mulch/exported-feature-table/table_P2BI_filtered.biom --table-type="OTU table" --to-hdf5 
   
   #import biom table back to qiime 2 artifact
       qiime tools import \
    --input-path table_P2BI_filtered.biom   \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path living_mulch_P2BI-filtered.qza
    
    #visualizing the resulting artifact 
       qiime feature-table summarize \
   --i-table living_mulch_P2BI-filtered.qza \
   --o-visualization living_mulch_P2BI-filtered.qzv
   
     #####create rarefraction table 
     qiime diversity alpha-rarefaction \
     --i-table living_mulch_P2BI-filtered.qza \
     --p-max-depth 30000 \
     --o-visualization 30000_rarefaction_table_for_living_mulch_P2BI-filtered.qzv
     
### according to the rarefraction curve, sequence depth of 20000 should capture most ASV. GA-20S-600 sample should be filtered out  

   #calculate alpha and beta diversity
     qiime diversity core-metrics-phylogenetic \
   --i-phylogeny rooted-tree.qza \
   --i-table living_mulch_P2BI-filtered.qza \
   --p-sampling-depth 20000 \
   --m-metadata-file living_mulch_data_metafile_massalin.tsv \
   --output-dir core-metrics-results-orginal-table-20000
   
############

### I tired to use decontam package in R to indentify possible contamination sequences caused by P2BS and P5BS blank 
#### however the resulting contanminated sequences are identified based on concentration of only two or three samples 
##### since P2BS and P5BS are closely resemble to other samples in beta diversity, exlude their seuqence will destroy the dataset 
##### therefore I decide to remove two blank in the downstream analysis without filtering their sequences from the whole datasets 
    
  
####################################################

   #after the filtration is done in pycharm, convert the filtered feature table from txt to biom blank sample and inrow sample was filtered 
  biom convert -i /home/hl46161/publish_living_mulch/exported-feature-table/table_no_blank_no_20_filtered.csv -o /home/hl46161/publish_living_mulch/exported-feature-table/table_no_blank_no_20_filtered.biom --table-type="OTU table" --to-hdf5 
   
     
   #import biom table back to qiime 2 artifact
          qiime tools import \
    --input-path table_no_blank_no_20_filtered.biom   \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path living_mulch_no_blank_no_20_filtered.qza

        #visualizing the resulting artifact 
                qiime feature-table summarize \
   --i-table living_mulch_no_blank_no_20_filtered.qza \
   --o-visualization living_mulch_no_blank_no_20_filtered.qzv
   
   ####classify the taxonomy of no blank no GA-20S-S600 datasets 
        qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification living_mulch_no_blank_filtered_taxonomy.qza
    
     ####filtered out the mitochondira and chloroplast sequences 
     qiime taxa filter-table \
  --i-table living_mulch_no_blank_no_20_filtered.qza \
  --i-taxonomy living_mulch_no_blank_filtered_taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast.qza
  
    ###### summarize the filtered ASV table This is the table I am going to use for downstream alpha and beta diversity analysis 
    qiime feature-table summarize \
   --i-table living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast.qza \
   --o-visualization living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast.qzv
   
  
   #calculate alpha and beta diversity
    qiime diversity core-metrics-phylogenetic \
   --i-phylogeny rooted-tree.qza \
   --i-table living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast.qza \
   --p-sampling-depth 19000 \
   --m-metadata-file living_mulch_data_metafile_massalin.tsv \
   --output-dir core-metrics-results-living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast-19000
  

#use a for loop to compare four alpha diversity metric between varaibles faith_pd_vector,evenness_vector,shannon_vector,observed_features_vector
# This loop need a document contain the name of alpha ddiversity metric  (alpha_diversity_metrics.txt)
#Also correlate alpha diversity metrics with the soil physical data 
#
while read line ; do
    echo $line
    path=core-metrics-results-living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast-19000
    
    qiime diversity alpha-group-significance \
  --i-alpha-diversity $path/${line}.qza \
  --m-metadata-file living_mulch_data_metafile_massalin.tsv \
  --o-visualization $path/${line}-group-sginificance.qzv
    
      qiime diversity alpha-correlation \
  --i-alpha-diversity $path/${line}.qza \
  --m-metadata-file living_mulch_data_metafile_massalin.tsv \
  --o-visualization $path/${line}_correlation.qzv
  
done <alpha_diversity_metrics.txt


#use a for loop to compare beta diversity 
## again this needs a document contain the beta diversity I would like to compare 
while read line ; do
    echo $line
    path=core-metrics-results-living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast-19000
    
    qiime diversity beta-group-significance \
  --i-distance-matrix $path/${line}.qza \
  --m-metadata-file living_mulch_data_metafile_massalin.tsv \
  --o-visualization $path/${line}_date_taken.qzv \
  --m-metadata-column Date_taken \
  --p-pairwise
 
   qiime diversity beta-group-significance \
  --i-distance-matrix $path/${line}.qza \
  --m-metadata-file living_mulch_data_metafile_massalin.tsv \
  --o-visualization $path/${line}_treatment.qzv \
  --m-metadata-column Treatment \
  --p-pairwise
  
  
done <beta_diversity_metrics.txt


##################### now we need to seperate the data into three sampling date and visualize the table 

         qiime feature-table filter-samples \
  --i-table living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast.qza \
  --m-metadata-file living_mulch_data_metafile_massalin.tsv \
  --p-where "[Date_taken] IN ('2018-05-21')" \
  --o-filtered-table May-table-no-mitochondria-no-chloroplast.qza
  
        qiime feature-table summarize \
   --i-table May-table-no-mitochondria-no-chloroplast.qza \
   --o-visualization May-table-no-mitochondria-no-chloroplast.qzv
  
           qiime feature-table filter-samples \
  --i-table living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast.qza \
  --m-metadata-file living_mulch_data_metafile_massalin.tsv \
  --p-where "[Date_taken] IN ('2018-06-28')" \
  --o-filtered-table June-table-no-mitochondria-no-chloroplast.qza
  
        qiime feature-table summarize \
   --i-table June-table-no-mitochondria-no-chloroplast.qza \
   --o-visualization June-table-no-mitochondria-no-chloroplast.qzv
  
        qiime feature-table filter-samples \
  --i-table living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast.qza \
  --m-metadata-file living_mulch_data_metafile_massalin.tsv \
  --p-where "[Date_taken] IN ('2018-08-31')" \
  --o-filtered-table August-table-no-mitochondria-no-chloroplast.qza
        
        qiime feature-table summarize \
   --i-table August-table-no-mitochondria-no-chloroplast.qza \
   --o-visualization August-table-no-mitochondria-no-chloroplast.qzv
  
##################################### do alpha diversity significance analysis to each sampling date 
### the sampling depth is determined by lowest sample depth in each of these subset data
    
    qiime diversity core-metrics-phylogenetic \
   --i-phylogeny rooted-tree.qza \
   --i-table May-table-no-mitochondria-no-chloroplast.qza \
   --p-sampling-depth 19000 \
   --m-metadata-file living_mulch_data_metafile_massalin.tsv \
   --output-dir core-metrics-results-May-table-no-mitochondria-no-chloroplast_filtered-19000
    
    qiime diversity core-metrics-phylogenetic \
   --i-phylogeny rooted-tree.qza \
   --i-table June-table-no-mitochondria-no-chloroplast.qza \
   --p-sampling-depth 28000 \
   --m-metadata-file living_mulch_data_metafile_massalin.tsv \
   --output-dir core-metrics-results-June-table-no-mitochondria-no-chloroplast_filtered-28000
   
       qiime diversity core-metrics-phylogenetic \
   --i-phylogeny rooted-tree.qza \
   --i-table August-table-no-mitochondria-no-chloroplast.qza \
   --p-sampling-depth 26000 \
   --m-metadata-file living_mulch_data_metafile_massalin.tsv \
   --output-dir core-metrics-results-August-table-no-mitochondria-no-chloroplast_filtered-26000

##################################### do alpha correplation analysis to each sampling date 

results_path="core-metrics-results-May-table-no-mitochondria-no-chloroplast_filtered-19000 core-metrics-results-June-table-no-mitochondria-no-chloroplast_filtered-28000 core-metrics-results-August-table-no-mitochondria-no-chloroplast_filtered-26000"
    
while read line ; do
    echo $line
    
    for path in $results_path
do
      echo $path
    
    qiime diversity alpha-group-significance \
  --i-alpha-diversity $path/${line}.qza \
  --m-metadata-file living_mulch_data_metafile_massalin.tsv \
  --o-visualization $path/${line}-group-sginificance.qzv
    
      qiime diversity alpha-correlation \
  --i-alpha-diversity $path/${line}.qza \
  --m-metadata-file living_mulch_data_metafile_massalin.tsv \
  --p-method spearman \
  --o-visualization $path/${line}_spearman_correlation.qzv
  
   qiime diversity alpha-correlation \
  --i-alpha-diversity $path/${line}.qza \
  --m-metadata-file living_mulch_data_metafile_massalin.tsv \
  --p-method pearson \
  --o-visualization $path/${line}_pearson_correlation.qzv
  
done
  
done <alpha_diversity_metrics.txt


################################### compared beta diveristy between treatment in each sampling date 

while read line ; do
    echo $line
    
      for path in $results_path
do
      echo $path
 
   qiime diversity beta-group-significance \
  --i-distance-matrix $path/${line}.qza \
  --m-metadata-file living_mulch_data_metafile_massalin.tsv \
  --o-visualization $path/${line}_treatment.qzv \
  --m-metadata-column treatment \
  --p-pairwise

done
  
done <beta_diversity_metrics.txt

##################################################### 

###beta diversity correlation with soil physical characteristic using mantel tests using both pearson correlation and spearman correlation
### each sampling date was done seperately 

results_path="core-metrics-results-May-table-no-mitochondria-no-chloroplast_filtered-19000 core-metrics-results-June-table-no-mitochondria-no-chloroplast_filtered-28000 core-metrics-results-August-table-no-mitochondria-no-chloroplast_filtered-26000"

correlation_method=("spearman pearson")

     for file in $results_path; \
do \
  echo $file
  
  
  while read line ; do
    echo $line
    
       qiime metadata distance-matrix \
  --m-metadata-file living_mulch_data_metafile_massalin_no_blank.tsv \
  --m-metadata-column $line \
  --o-distance-matrix $file/${line}_distance_matrix.qza 
  
for method in $correlation_method; \
  do 
      echo $method

  qiime diversity mantel \
  --i-dm1 $file/{$line}_distance_matrix.qza \
  --i-dm2 $file/weighted_unifrac_distance_matrix.qza \
  --p-method $method \
  --p-permutations 999 \
  --p-label1 ${line}_concentration \
  --p-label2 weighted_unifrac_distance \
  --p-intersect-ids TRUE \
  --o-visualization $file/${line}_weighted_unifrac_spearman_association.qzv
  
    qiime diversity mantel \
  --i-dm1 $file/{$line}_distance_matrix.qza  \
  --i-dm2 $file/unweighted_unifrac_distance_matrix.qza \
  --p-method $method \
  --p-permutations 999 \
  --p-label1 $line \
  --p-label2 unweighted_unifrac_distance \
  --p-intersect-ids TRUE \
  --o-visualization $file/${line}_unweighted_unifrac_spearman_association.qzv
  
      qiime diversity mantel \
  --i-dm1 $file/{$line}_distance_matrix.qza  \
  --i-dm2 $file/bray_curtis_distance_matrix.qza \
  --p-method $method \
  --p-permutations 999 \
  --p-label1 $line \
  --p-label2 bray_curtis_distance \
  --p-intersect-ids TRUE \
  --o-visualization $file/{$line}_bray_curtis_spearman_association.qzv
  
done
  
done <column_name_list.txt

done

  
#####################################################################
  
  #export the no mitochondira and no chloroplast table out for differential abudance test in deseq2 
  
   qiime tools export \
    --input-path living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast.qza \
    --output-path exported-feature-table
    
   biom convert -i feature-table.biom -o living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast.txt --to-tsv
   
##########################################################################

  #export the no mitochondira and no chloroplast table out for NMDS analysis in R 
  
   qiime tools export \
    --input-path /home/hl46161/publish_living_mulch/core-metrics-results-living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast-19000/rarefied_table.qza \
    --output-path /home/hl46161/publish_living_mulch/exported-feature-table-for-NMDS/
    
    cd /home/hl46161/publish_living_mulch/exported-feature-table-for-NMDS/
    
   biom convert -i feature-table.biom -o living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast-19000-rarefied.txt --to-tsv
   

   
##############################################################################

         qiime feature-table filter-samples \
  --i-table living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast.qza \
  --m-metadata-file living_mulch_data_metafile_massalin.tsv \
  --p-where "[treatment] IN ('CerealRye')" \
  --o-filtered-table CerealRye-table-no-mitochondria-no-chloroplast.qza
  
    qiime feature-table summarize \
   --i-table CerealRye-table-no-mitochondria-no-chloroplast.qza \
   --o-visualization CerealRye-table-no-mitochondria-no-chloroplast.qzv
  
         qiime feature-table filter-samples \
  --i-table living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast.qza \
  --m-metadata-file living_mulch_data_metafile_massalin.tsv \
  --p-where "[treatment] IN ('CrimsonClover')" \
  --o-filtered-table CrimsonClover-table-no-mitochondria-no-chloroplast.qza
  
    qiime feature-table summarize \
   --i-table CrimsonClover-table-no-mitochondria-no-chloroplast.qza \
   --o-visualization CrimsonClover-table-no-mitochondria-no-chloroplast.qzv
  
         qiime feature-table filter-samples \
  --i-table living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast.qza \
  --m-metadata-file living_mulch_data_metafile_massalin.tsv \
  --p-where "[treatment] IN ('LivingMulch')" \
  --o-filtered-table LivingMulch-table-no-mitochondria-no-chloroplast.qza
  
      qiime feature-table summarize \
   --i-table LivingMulch-table-no-mitochondria-no-chloroplast.qza \
   --o-visualization LivingMulch-table-no-mitochondria-no-chloroplast.qzv
  
  
           qiime feature-table filter-samples \
  --i-table living_mulch_no_blank_no_20_filtered-no-mitochondria-no-chloroplast.qza \
  --m-metadata-file living_mulch_data_metafile_massalin.tsv \
  --p-where "[treatment] IN ('NoCover')" \
  --o-filtered-table NoCover-table-no-mitochondria-no-chloroplast.qza
  
      qiime feature-table summarize \
   --i-table NoCover-table-no-mitochondria-no-chloroplast.qza \
   --o-visualization NoCover-table-no-mitochondria-no-chloroplast.qzv

##############################################################################

   qiime diversity core-metrics-phylogenetic \
   --i-phylogeny rooted-tree.qza \
   --i-table CerealRye-table-no-mitochondria-no-chloroplast.qza \
   --p-sampling-depth 25000 \
   --m-metadata-file living_mulch_data_metafile_massalin.tsv \
   --output-dir core-metrics-results-CerealRye-table-no-mitochondria-no-chloroplast-25000
   
    qiime diversity core-metrics-phylogenetic \
   --i-phylogeny rooted-tree.qza \
   --i-table CrimsonClover-table-no-mitochondria-no-chloroplast.qza \
   --p-sampling-depth 19000 \
   --m-metadata-file living_mulch_data_metafile_massalin.tsv \
   --output-dir core-metrics-results-CrimsonClover-table-no-mitochondria-no-chloroplast-19000
   
   
    qiime diversity core-metrics-phylogenetic \
   --i-phylogeny rooted-tree.qza \
   --i-table  LivingMulch-table-no-mitochondria-no-chloroplast.qza \
   --p-sampling-depth 31000 \
   --m-metadata-file living_mulch_data_metafile_massalin.tsv \
   --output-dir core-metrics-results-LivingMulch-table-no-mitochondria-no-chloroplast-31000
   
    qiime diversity core-metrics-phylogenetic \
   --i-phylogeny rooted-tree.qza \
   --i-table NoCover-table-no-mitochondria-no-chloroplast.qza \
   --p-sampling-depth 34000 \
   --m-metadata-file living_mulch_data_metafile_massalin.tsv \
   --output-dir core-metrics-results-NoCover-table-no-mitochondria-no-chloroplast-34000

###################################################################################################
  
   results_path="core-metrics-results-NoCover-table-no-mitochondria-no-chloroplast-34000 core-metrics-results-LivingMulch-table-no-mitochondria-no-chloroplast-31000 core-metrics-results-CrimsonClover-table-no-mitochondria-no-chloroplast-19000 core-metrics-results-CerealRye-table-no-mitochondria-no-chloroplast-25000"
    
while read line ; do
    echo $line
    
    for path in $results_path
do
      echo $path
    
    qiime diversity alpha-group-significance \
  --i-alpha-diversity $path/${line}.qza \
  --m-metadata-file living_mulch_data_metafile_massalin.tsv \
  --o-visualization $path/${line}-group-sginificance.qzv
    
      qiime diversity alpha-correlation \
  --i-alpha-diversity $path/${line}.qza \
  --m-metadata-file living_mulch_data_metafile_massalin.tsv \
  --o-visualization $path/${line}_correlation.qzv
done
  
done <alpha_diversity_metrics.txt


################################################




