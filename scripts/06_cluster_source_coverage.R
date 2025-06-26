
library(GenomicRanges)
library(dplyr)
library(openxlsx)
library(Biostrings)
setwd("~/integration_site/results")
all_mice <- list.files(pattern="^\\d+\\.bycoords\\.bam$")
all_mice <-  gsub("\\D","", all_mice) |> as.numeric()


#d <- read.xlsx("isa.audit.jd862.xlsx") |> filter(classification) |> filter(!is.na(source))
#d$source_coverage <- NA



sources <- c('chr3:54756601-54763754','chr3:52137362-52140547','chr13:100725460-100730733','chr3:60397294-60402505','chr11:116280231-116285535','chr4:126902313-126909484','chr11:116676420-116683530','chr18:47781821-47787021','chr2:7408884-7414163','chrX:22430173-22437344','chr1:142512630-142517888','chr6:31776305-31781572','chr4:140934913-140937903','chr11:105919072-105924386','chr16:31966895-31972157','chr3:156505180-156510398','chr2:72100706-72105943','chr17:45955573-45960846','chr3:139825771-139831055','chr6:105007550-105010276','chr10:96028398-96035564','chr2:118385247-118388861','chr7:82915270-82922443','chr18:68545315-68552467','chr10:116442165-116449281','chr8:45895952-45899881','chr13:50885311-50889705','chr1:28495216-28498223','chr4:40054170-40061404','chr18:27131386-27138509','chr6:28924174-28929352','chr1:108336585-108343768','chr7:4458441-4463701','chr17:54122694-54129891','chr10:77786566-77793719','chrX:90860584-90867717','chr3:82388101-82395311','chr14:64700829-64706164','chr12:32078022-32082383','chr10:59217621-59222882','chr8:117729304-117734453','chr5:81525547-81532649','chr5:99959513-99962262','chr13:85747029-85753599','chr10:78652305-78658886','chr6:38253374-38259959','chr1:81541884-81548964','chr2:154196183-154201512')


get_source_coverage <- function(virus, Expt.ID) {
  regmatches(virus, regexec(r"(^(chr\w+):(\d+)-(\d+)$)",virus))[[1]] -> virus.re
  virus.gr <- GRanges(virus.re[2],IRanges(start=as.numeric(virus.re[3]),end=as.numeric(virus.re[4])))
  if (width(virus.gr)>10000) return(0)
  r <- strsplit(system(paste0("/rds/project/sjl83/rds-sjl83-green/jd862/mambaforge/envs/samtools/bin/samtools view ",Expt.ID,".bycoords.bam ",virus), intern=T),"\t")


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
if (file.exists("source.coverage.xlsx")) {
  d <- read.xlsx("source.coverage.xlsx")
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

  write.xlsx(d,"source.coverage.xlsx")
}

d <- read.xlsx("isa.audit.jd862.xlsx")
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
 as.data.frame(d) |> openxlsx::write.xlsx("isa.audit.jd862.source_coverage.xlsx")
}


