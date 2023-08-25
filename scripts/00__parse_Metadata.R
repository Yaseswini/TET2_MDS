## Metadata for the TET2 MDS project ## 

options( stringsAsFactors = F )
set.seed(41568)

library( tidyverse )
library( rtracklayer )
library( plotrix )
library( magrittr )
library( methylKit )
library( ggthemes )
library( ggrepel )
library( grid )
library( gridExtra )
library( VennDiagram )
library( xlsx )

source( "./scripts/functions/annotate_loci_old.R" )
source( "./scripts/functions/05_RunEnrichR.R" )
source( "./scripts/functions/makeDataFramefromGRanges.R" )

## === Metadata == ## 
### Read methylcall files from GEO submission
### Download the methylcall files from GEO accession ID : 

methylcallFileList = list.files( <GEO_dir> , full.names = T , pattern = "*_methylcall.txt" )
names(methylcallFileList) = gsub( "_methylcall.txt","",basename(methylcallFileList) )

ERRBS_metadata = data.frame( methylcallFile = methylcallFileList , methylcallFilename = names(methylcallFileList) ) %>% 
                            remove_rownames() %>% 
                            rowwise() %>%
                            mutate( group = case_when( gsub("ERRBS_|_rep[1|2]","",methylcallFilename) == "15month_tet2_mut" ~ "TET15" ,
                                                       gsub("ERRBS_|_rep[1|2]","",methylcallFilename) == "15month_tet2_wt" ~ "WT15" ,
                                                       gsub("ERRBS_|_rep[1|2]","",methylcallFilename) == "4month_tet2_mut" ~ "TET4" ,
                                                       gsub("ERRBS_|_rep[1|2]","",methylcallFilename) == "4month_tet2_wt" ~ "WT4" ) )

### Mapping gene ids from the rnaseq data to their TSS sites 
rnaseq_counts %>% dplyr::select( Chr , TSS , TSE , ID , Gene )
RNAseq_gene_annotationData = read_tsv( "./data/referenceGenome/" ) 
RNAseq_gene_annotationData = RNAseq_gene_annotationData %>% mutate( geneLength = TSE - TSS + 1 ) %>%
                                      dplyr::select( -Chr , -TSS , -TSE , -Gene ) %>% 
                                      plyr::rename( c('ID'='gene_id') ) %>% 
                                      remove_rownames() %>%
                                      column_to_rownames( var = "gene_id" )

## ======= Genome annotation ===== ## 
### Reading GTF file 
Zv9_gtfData = readGFF( "./data/referenceGenome/gtfFile" )
Zv9_gtfSimple = Zv9_gtfData %>% 
                        dplyr::select( seqid , gene_id , gene_name , strand , gene_biotype ) %>% 
                        unique %>% 
                        dplyr::rename( "Zebrafish_ensid" = "gene_id" , "Zebrafish_symbol" = "gene_name" )

### Reading the orthology data 
Zv9_orthologData = read.table( orthologyFile ,head = T , sep = "\t" ) %>% 
                              filter( Rank == 1 | Rank == 2 ) %>% 
                              dplyr::rename( "Zebrafish_ensid" = "Zebrafish_ID" ) %>% 
                              separate_rows( Other_possibilities )
Zv9_orthologData = rbind( setNames( Zv9_orthologData[ , c("Zebrafish_ensid","Orthologous_human_gene") ] , c("Zebrafish_ensid","HumanGene") ) , 
                          setNames( Zv9_orthologData[ , c("Zebrafish_ensid","Other_possibilities") ] , c("Zebrafish_ensid","HumanGene") ) )
Zv9_orthologData = Zv9_orthologData %>% filter( HumanGene != "" )
Zv9_orthologData = unique( Zv9_orthologData )
Zv9_orthologData = Zv9_orthologData  %>% group_by( Zebrafish_ensid ) %>% dplyr::summarize( HumanGene = paste0( unique(HumanGene) , collapse=";" ) )

Zv9_proteinCoding_ENSIDs = Zv9_gtfSimple %>% filter( gene_biotype == "protein_coding" ) 

geneAnnot.rnaseq = rnaseq_counts %>% 
                        dplyr::select( Chr , TSS , TSE , ID ) %>% 
                        mutate( symbol = gsub("!.*","",ID) ) %>% 
                        dplyr::rename( "Zebrafish_ensid" = "symbol" , 
                                      "geneID" = "ID" , 
                                      "chr" = "Chr" , 
                                      "gene_start" = "TSS" , 
                                      "gene_end" = "TSE" ) %>%
                        left_join( Zv9_gtfSimple %>% 
                                        dplyr::select( Zebrafish_ensid , Zebrafish_symbol , strand , gene_biotype ) ) %>%
                        mutate( TSS = if_else( strand == "+" , gene_start , gene_end ) ) %>% 
                        mutate( promoter_start = TSS - 500 , promoter_end = TSS + 500 ) %>% 
                        left_join( Zv9_orthologData )

promLoci_500_and_250_bp = rnaseq_counts %>% 
                            dplyr::select( Chr , TSS , TSE , ID ) %>% 
                            mutate( symbol = gsub("!.*","",ID) ) %>% 
                            dplyr::rename( "Zebrafish_ensid" = "symbol" , 
                                           "geneID" = "ID" , 
                                           "chr" = "Chr" , 
                                           "gene_start" = "TSS" , 
                                           "gene_end" = "TSE" ) %>%
                            left_join( Zv9_gtfSimple %>% 
                                          dplyr::select( Zebrafish_ensid , Zebrafish_symbol , strand , gene_biotype ) ) %>%
                            mutate( TSS = if_else( strand == "+" , gene_start , gene_end ) ) %>% 
                            mutate( promoter_start = TSS - 500 , promoter_end = TSS + 250 ) %>% 
                            left_join( Zv9_orthologData ) %>% 
                            dplyr::select( chr , promoter_start , promoter_end , strand , Zebrafish_ensid , Zebrafish_symbol , HumanGene , gene_biotype ) %>% 
                            dplyr::rename( "start" = "promoter_start" , 
                                           "end" = "promoter_end" , 
                                           "id" = "Zebrafish_ensid" , "symbol" = "Zebrafish_symbol" ) %>%
                            unique 

#promLoci_500_and_250_bp = promLoci_500_and_250_bp %>% filter( gene_biotype == "protein_coding" )
promLoci_500_and_250_bp_gr = makeGRangesFromDataFrame( unique(promLoci_500_and_250_bp)  , keep.extra.columns = T )
promoterDF = makeDataFramefromGRanges( promLoci_500_and_250_bp_gr ) 
promoterDF = promoterDF %>% mutate( "promoter_loci" = paste0( chr , ":" , start , "!" , end ) ) %>% mutate( geneID = paste0( id , "!" , symbol) )



## CpGislands 
CpGisland_gr = readRDS( "./data/referenceGenome/Zv9_CpGi.rds" )
CpGisland_df = makeDataFramefromGRanges( CpGisland_gr ) %>% 
                   mutate( chrBase = paste0( chr , "." , start , "." , end ) ) %>% 
                    dplyr::rename( "CpG_id" = "symbol" ) %>% 
                    dplyr::select( -id , -strand )
Map_CpGi_to_genomicAnnot = annotate_sites( loci = CpGisland_df[, c("chr","start","end") ] , 
                               referenceAnnot = genomeLoci , 
                               preference = c("Promoter","Exon","firstIntron","Intron","Intergenic") , 
                               overlap.type = "within" , 
                               min.overlap = 100L ) %>% 
                          merge_annotationRes() %>%
                          left_join( CpGisland_df %>% dplyr::select( chrBase , CpG_id ) , by = "chrBase" )
Map_CpGi_to_genomicAnnot = Map_CpGi_to_genomicAnnot %>% dplyr::select( chrBase , CpG_id , id , symbol , Region ) %>% unique %>% dplyr::rename( "CpGi_loci" = "chrBase" )


CpGisland_mapped_withinPromoter = annotate_sites( loci = CpGisland_df[, c("chr","start","end") ] , 
                                                   referenceAnnot = list( "Promoter" = makeGRangesFromDataFrame( promoterDF , keep.extra.columns = T ) ) , 
                                                   preference = "Promoter" , 
                                                   overlap.type = "within" , 
                                                   min.overlap = 100L ) %>% 
                                                  merge_annotationRes() %>%
                                                  left_join( CpGisland_df %>% dplyr::select( chrBase , CpG_id ) , by = "chrBase" ) %>% unique %>% 
                                                  filter( Region == "Promoter" )
CpGisland_mapped_withinPromoter = CpGisland_mapped_withinPromoter %>% dplyr::select( chrBase , CpG_id , id , symbol , Region ) %>% unique %>% dplyr::rename( "CpGi_loci" = "chrBase" )


