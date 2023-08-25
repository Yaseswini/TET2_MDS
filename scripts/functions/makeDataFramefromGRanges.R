makeDataFramefromGRanges <- function( grangesObj )
{
  dat = data.frame( chr = seqnames(grangesObj) ,
                    start = start( grangesObj ) ,
                    end = end( grangesObj ) ,
                    strand = strand( grangesObj ) ,
                    elementMetadata(grangesObj) )
  return(dat)
}