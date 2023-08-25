## Perform differential methylation ## 

comparisons = c("TET4.vs.WT4","TET15.vs.TET4","WT15.vs.WT4","TET15.vs.WT15")

methylkit_dmr_analysis = list()
methDiff_thresh = 15 ; 
qvalue_thresh = 0.05 
colors = c( "Promoter" = "#20854EFF" , "Hyper" = "#DF8F44FF" , "Hypo" = "#0073C2FF" , "Unchanged" = "grey" )

for( idx in 1:length(comparisons) ){

	comp.split = unlist(strsplit(comparisons[[idx]],".vs."))
	comp1.rep = ERRBS_metadata %>% filter( group == comp.split[[1]] ) %>% .$methylcallFile 
	comp1.ids = ERRBS_metadata %>% filter( group == comp.split[[1]] ) %>% .$methylcallFilename
	comp2.rep = ERRBS_metadata %>% filter( group == comp.split[[2]] ) %>% .$methylcallFile 
	comp2.ids = ERRBS_metadata %>% filter( group == comp.split[[2]] ) %>% .$methylcallFilename


	sample.id = c(comp1.ids,basename(comp2.ids))
	treatment = c(rep(1,length(comp1.rep)),rep(0,length(comp2.rep)))
	assembly = "Zv9"
	fileList = c( comp1.rep , comp2.rep )
	context = "CpG"

	myobj = methRead( as.list(fileList) , sample.id = as.list( sample.id ) , assembly=assembly , treatment=treatment , context=context)
		
		## ----- Promoter ----- ## 
		myobj_promoter = regionCounts( myobj , promLoci_500_and_250_bp_gr )
		myobj_promoter_filt <- filterByCoverage(myobj_promoter,
						                      lo.count=10,
						                      lo.perc=NULL,
						                      hi.count=NULL,
						                      hi.perc=90)

		# Perform simple normalization
		myobj_promoter_filt_norm <- normalizeCoverage(myobj_promoter_filt, method = "mean")
		# Merge the samples again
		meth_promoter <- unite(myobj_promoter_filt_norm)

		# Test for differential methylation... This might take a few minutes.
		myDiff_promoter <- calculateDiffMeth(meth_promoter,adjust="BH")

		myDiff_promoter_dataFrame = data.frame( myDiff_promoter ) %>% mutate( promoter_loci = paste0( chr , ":" , start , "!" , end ) )
		myDiff_promoter_dataFrame = myDiff_promoter_dataFrame %>% unique %>%
												mutate( methGroup = case_when( meth.diff > methDiff_thresh & qvalue < qvalue_thresh ~ "Hyper" , 
																			   meth.diff < -methDiff_thresh & qvalue < qvalue_thresh ~ "Hypo" ,
																			   qvalue > 0 | abs(meth.diff) < methDiff_thresh & qvalue < 0 ~ "NS" ) )
		myDiff_promoter_dataFrame = myDiff_promoter_dataFrame %>% left_join( promoterDF %>% dplyr::select( -chr,-start,-end,-strand) , by = "promoter_loci" ) %>% unique

		myDiff_promoter_dataFrame_tally = myDiff_promoter_dataFrame %>% 
												dplyr::select( chr,start,end,pvalue,qvalue,meth.diff ) %>% 
												unique %>%
												mutate( methGroup = case_when( meth.diff > methDiff_thresh & qvalue < qvalue_thresh ~ "Hyper" , 
																			   meth.diff < -methDiff_thresh & qvalue < qvalue_thresh ~ "Hypo" ,
																			   qvalue > 0 | ( abs(meth.diff) < methDiff_thresh & qvalue < 0 ) ~ "NS" ) ) %>%
												group_by( methGroup ) %>%
												tally()

		hypermethyl = myDiff_promoter_dataFrame %>% filter( methGroup == "Hyper" ) %>% unique 
		hypomethyl =  myDiff_promoter_dataFrame %>% filter( methGroup == "Hypo" ) %>% unique
		length( unique(hypermethyl[["promoter_loci"]]) )
		length( unique(hypomethyl[["promoter_loci"]]) )

		hypermethyl_humangenes = sort( unique( as.vector( na.omit( unlist( strsplit( hypermethyl[["HumanGene"]] , ";" ) ) ) ) ) )
		hypomethyl_humangenes = sort( unique( as.vector( na.omit( unlist( strsplit( hypomethyl[["HumanGene"]] , ";" ) ) ) ) ) )

		enrichResult_hypermethyl = run_enrichR_human( hypermethyl_humangenes , outDIR = comp_dir, outLabel = paste0( comparisons[[idx]] , "_HyperMethylatedPromoters" ) ) 
		enrichResult_hypermethyl_sig = enrichResult_hypermethyl %>% filter( P.value < 0.05 )
		enrichResult_hypermethyl_sigList = split( enrichResult_hypermethyl_sig  , enrichResult_hypermethyl_sig$category )
     	enrichResult_hypermethyl_sigList = lapply( enrichResult_hypermethyl_sigList , function(df) df %>% arrange( P.value ) )

		Sys.sleep(60)

		enrichResult_hypomethyl = run_enrichR_human( hypomethyl_humangenes , outDIR = comp_dir , outLabel = paste0( comparisons[[idx]], "_HypoMethylatedPromoters" ) ) 
		enrichResult_hypomethyl_sig = enrichResult_hypomethyl %>% filter( P.value < 0.05 )
		enrichResult_hypomethyl_sigList = split( enrichResult_hypomethyl_sig  , enrichResult_hypomethyl_sig$category )
     	enrichResult_hypomethyl_sigList = lapply( enrichResult_hypomethyl_sigList , function(df) df %>% arrange( P.value ) )

     	Sys.sleep(60)

		methylkit_dmr_analysis[[ comparisons[[idx]] ]][[ "myDiff_promoter" ]] = myDiff_promoter_dataFrame
		methylkit_dmr_analysis[[ comparisons[[idx]] ]][[ "myDiff_promoter_tally" ]] = myDiff_promoter_dataFrame_tally
		methylkit_dmr_analysis[[ comparisons[[idx]] ]][[ "HyperDMRs_enrichR_unfiltered" ]] = enrichResult_hypermethyl
		methylkit_dmr_analysis[[ comparisons[[idx]] ]][[ "HypoDMRs_enrichR_unfiltered" ]] = enrichResult_hypomethyl
		methylkit_dmr_analysis[[ comparisons[[idx]] ]][[ "HyperDMRs_enrichR" ]] = enrichResult_hypermethyl_sigList
		methylkit_dmr_analysis[[ comparisons[[idx]] ]][[ "HypoDMRs_enrichR" ]] = enrichResult_hypomethyl_sigList
		methylkit_dmr_analysis[[ comparisons[[idx]] ]][[ "HyperDMRs_MDSgeneset_overlap" ]] = lapply( MDS_gmtData , function(x) intersect( x , hypermethyl_humangenes ) )
		methylkit_dmr_analysis[[ comparisons[[idx]] ]][[ "HypoDMRs_MDSgeneset_overlap" ]] = lapply( MDS_gmtData , function(x) intersect( x , hypomethyl_humangenes ) )
		
		
		## ----- CpGislands ----- ## 
		CpGisland_df_1 = CpGisland_df %>% mutate( "CpGi_loci" = paste0( chr,".",start,".",end) )
		myobj_CpGislands = regionCounts( myobj , CpGisland_gr  )
		myobj_CpGislands_filt <- filterByCoverage( myobj_CpGislands ,
						                      		lo.count=10,
						                      		lo.perc=NULL,
						                      		hi.count=NULL,
						                      		hi.perc=90 )

		# Perform simple normalization
		myobj_CpGislands_filt_norm <- normalizeCoverage( myobj_CpGislands_filt , method = "mean" )
		# Merge the samples again
		meth_CpGislands <- unite(myobj_CpGislands_filt_norm)

		# Test for differential methylation... This might take a few minutes.
		myDiff_CpGislands <- calculateDiffMeth( meth_CpGislands , adjust="BH" )

		myDiff_CpGislands_dataFrame = data.frame( myDiff_CpGislands ) %>% mutate( CpGi_loci = paste0( chr , "." , start , "." , end ) )
		myDiff_CpGislands_dataFrame = myDiff_CpGislands_dataFrame %>% unique %>%
												mutate( methGroup = case_when( meth.diff > methDiff_thresh & qvalue < qvalue_thresh ~ "Hyper" , 
																			   meth.diff < -methDiff_thresh & qvalue < qvalue_thresh ~ "Hypo" ,
																			   qvalue > 0 | ( abs(meth.diff) < methDiff_thresh & qvalue < 0 ) ~ "NS" ) )
		myDiff_CpGislands_dataFrame = myDiff_CpGislands_dataFrame %>% left_join( CpGisland_df_1 %>% dplyr::select( CpGi_loci , CpG_id ) , by = "CpGi_loci" ) %>% unique

		myDiff_CpGislands_dataFrame_tally = myDiff_CpGislands_dataFrame %>% 
												dplyr::select( chr,start,end,pvalue,qvalue,meth.diff ) %>% 
												unique %>%
												mutate( methGroup = case_when( meth.diff > methDiff_thresh & qvalue < qvalue_thresh ~ "Hyper" , 
																			   meth.diff < -methDiff_thresh & qvalue < qvalue_thresh ~ "Hypo" ,
																			   qvalue > 0 | abs(meth.diff) < methDiff_thresh & qvalue < 0 ~ "NS" ) ) %>%
												group_by( methGroup ) %>%
												tally()

		## Annotating CpGislands to genomic regions 
		myDiff_CpGislands_dataFrame_annotToGene = myDiff_CpGislands_dataFrame %>% left_join( Map_CpGi_to_genomicAnnot , by = c( "CpGi_loci" , "CpG_id" ) )

		## === Combine all results === ## 
		methylkit_dmr_analysis[[ comparisons[[idx]] ]][[ "myDiff_CpGislands" ]] = myDiff_CpGislands_dataFrame
		methylkit_dmr_analysis[[ comparisons[[idx]] ]][[ "myDiff_CpGislands_tally" ]] = myDiff_CpGislands_dataFrame_tally
		methylkit_dmr_analysis[[ comparisons[[idx]] ]][[ "CpGislands_annotToGenes" ]] = myDiff_CpGislands_dataFrame_annotToGene

		write_tsv( myDiff_promoter_dataFrame , file = file.path( comp_dir , paste0( comparisons[[idx]] , "_PromoterMethylationDifferences.tsv" ) ) )
		write_tsv( myDiff_CpGislands_dataFrame , file = file.path( comp_dir , paste0( comparisons[[idx]] , "_CpGislandsMethylationDifferences.tsv" ) ) )
}
