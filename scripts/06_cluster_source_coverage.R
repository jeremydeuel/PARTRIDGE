
library(GenomicRanges)
library(dplyr)
library(openxlsx)
library(Biostrings)
all_mice <- list.files("test_data",pattern="^\\w+\\.bam$")
all_mice <-  gsub("\\.bam","", all_mice)


#d <- read.xlsx("isa.audit.jd862.xlsx") |> filter(classification) |> filter(!is.na(source))
#d$source_coverage <- NA



sources <- c('chr6:31776305-31781572')


get_source_coverage <- function(virus, Expt.ID) {
  regmatches(virus, regexec(r"(^(chr\w+):(\d+)-(\d+)$)",virus))[[1]] -> virus.re
  virus.gr <- GRanges(virus.re[2],IRanges(start=as.numeric(virus.re[3]),end=as.numeric(virus.re[4])))
  if (width(virus.gr)>10000) return(0)
  r <- strsplit(system(paste0("samtools view test_data/",Expt.ID,".bam ",virus), intern=T),"\t")


  ref_cigar_len <- sapply(sapply(strsplit(sapply(r,`[[`,6), "(?<=[^0123456789])",perl=T),
                                 function(x) {
                                   return(sapply(x, function(y) {
                                     if (length(y)==0) return(0)
                                     if (substr(y,nchar(y),nchar(y))%in%c("M","D","N","=","X")) {
                                       return(as.numeric(substr(y,1,nchar(y)-1)))
                                     } else {
                                       return(0)
                                     }
                                   }))}), sum)

  synth_barcode <- nchar(sapply(r,`[[`,1))<46
  final_barcodes <- rep(NA,length(r))
  real_barcodes <- sapply(r,`[[`,1)[!synth_barcode]
  final_barcodes[!synth_barcode] <- substr(real_barcodes,nchar(real_barcodes)-8,nchar(real_barcodes))
  rm(real_barcodes)
  barcodes <- tibble(
    rname=sapply(r,`[[`,1)[synth_barcode],
    seq=sapply(r, `[[`, 10)[synth_barcode],
    rev=as.numeric(sapply(r, `[[`, 2))[synth_barcode],
  )
  if (nrow(barcodes)>0) {
    print("synthetising barcodes")
    barcodes$first <- bitwAnd(barcodes$rev,0x40)>0
    barcodes$rev <- !xor(bitwAnd(barcodes$rev,0x40)>0,bitwAnd(barcodes$rev,0x10)>0)
    barcodes$seq[!barcodes$rev] <- substr(barcodes$seq[!barcodes$rev],1,4)
    barcodes$seq[barcodes$rev] <- as.character(reverseComplement(DNAStringSet(substr(barcodes$seq[barcodes$rev],nchar(barcodes$seq[barcodes$rev])-4,nchar(barcodes$seq[barcodes$rev])))))
    barcodes$barcode_fwd <- ""
    barcodes$barcode_rev <- ""
    barcodes$barcode_fwd[barcodes$first] <- barcodes$seq[barcodes$first]
    barcodes$barcode_rev[!barcodes$first] <- barcodes$seq[!barcodes$first]
    barcodes |> group_by(rname) |> summarise(barcode=paste0(barcode_fwd,barcode_rev,collapse="")) -> barcodes
    final_barcodes[synth_barcode] <- barcodes$barcode[match(sapply(r,`[[`,1)[synth_barcode],barcodes$rname)]
    print(head(barcodes))
    rm(barcodes)
  }
  rm(synth_barcode)


  GRanges(sapply(r,`[[`,3),IRanges(start=as.numeric(sapply(r,`[[`,4)), end= ref_cigar_len + as.numeric(sapply(r,`[[`,4)) ),barcode=final_barcodes) -> reads
  reads[as.numeric(sapply(r,`[[`,5))>39] -> reads #MAPQ Filtering

  left <- virus.gr-25
  right <- virus.gr+25
  start(right) <- end(right)-50
  end(left) <- start(left)+50
  virus.gr <- c(left, right)

  return(reads$barcode[subjectHits(findOverlaps(virus.gr, reads, type="within", select="all"))] |> unique() |> length())
}
if (file.exists("test_data/source.coverage.xlsx")) {
  d <- read.xlsx("test_data/source.coverage.xlsx")
} else {
  d <- NULL
}


for (s in sources) {
  if (s %in% d$source) next
  try({
    for (m in all_mice) {
      print(paste0(s, ", ", m))
      d <- rbind(d, c(source=s, Expt.ID=m, coverage=get_source_coverage(s, m)))
    }
    d <- as.data.frame(d)
  })

  write.xlsx(d,"test_data/source.coverage.xlsx")
}

d <- read.xlsx("test_data/isa.audit.xlsx")
d$source_coverage <- NA

for (i in seq_len(nrow(d))) {
  if (!d$classification[i]) next
  if (!is.na(d$source_coverage[i])) next
  print(d[i,"name"])
  virus <- d$source[i]

  if (d$mouse[i]=="multi") {
    for (Expt.ID in all_mice) {
      print(paste0("... mouse ",Expt.ID))
      row <- d[i,]
      row$mouse <- Expt.ID
      try({
        row$source_coverage <-get_source_coverage(virus, Expt.ID)
      })

      d <- rbind(d,row)
    }
    d$source_coverage[i] <- Inf
  } else {
    print(paste0("... mouse ",d$mouse[i]))
    try({
      d$source_coverage[i] <- get_source_coverage(virus, d$mouse[i])
    })
  }
 as.data.frame(d) |> openxlsx::write.xlsx("test_data/isa.audit.source_coverage.xlsx")
}


