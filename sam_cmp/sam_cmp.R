library(Rsamtools)
library(stringr)
library(ggplot2)
library(reshape2)


ask_cigar <- function(cigar, op) {
  flag = 0
  if (length(cigar) == 1) {
    cigar <- str_extract_all(cigar, "\\*|([0-9]+)|([MIDNSHPX=])", simplify = T)
  }
  for (i in seq(from=2, to=length(cigar), by=2)) {
    if (cigar[i] == op) {
      flag = 1
      break
    }
  }
  if (flag == 1) {
    cigar[i - 1]
  } else {
    FALSE
  }
}

sum_cigar <- function(cigar_splited) {
  sum <- 0
  for (i in seq(from=2, to=length(cigar_splited), by=2)) {
    if ((cigar_splited[i] != "") & (cigar_splited[i] != "S")){
      sum <- sum + as.integer(cigar_splited[i - 1])
    }
  }
  sum
}



sam_data <- scanBam("SampleName.fixed.sorted.uniq.rg.bam", param=ScanBamParam(what=c("pos", "cigar", "qwidth", "qual")))
vcf <- read.delim("results.table_Wed#May#20#19-06-18#2020_MRPAS50", header = TRUE, sep = "\t")
what=scanBamWhat()
param=ScanBamParam(what=scanBamWhat())

end <- length(sam_data[[1]]$cigar)   #5000
#cigar_splited <- str_extract_all(sam_data[[1]]$cigar, "\\*|([0-9]+[MIDNSHPX=])", simplify = T)
cigar_splited <- str_extract_all(sam_data[[1]]$cigar[1:end], "\\*|([0-9]+)|([MIDNSHPX=])", simplify = T)

lm_pos <- sam_data[[1]]$pos[1:end]
rm_pos <- sam_data[[1]]$pos[1:end] + apply(cigar_splited, 1, sum_cigar) - 1


#reported_pos <- vcf$POS.Data[!is.na(vcf$POS.Data)]
reported_pos <- vcf$POS.Data[which(vcf$Positive == TRUE)]
idxs <- lapply(reported_pos, function(pos, l_pos, r_pos) {which(l_pos <= pos & pos <= r_pos)}, l_pos = lm_pos, r_pos = rm_pos)

pp <- c()
l <- c()
r <- c()

p <- c()
q <- c()


for (i in 1 : length(idxs)) {
  l <- length(idxs[[i]])
  if (l > 0) {
    pos <- integer(length = l)
    qual <- integer(length = l)
    for (j in 1 : l) {
      pos[j] <- reported_pos[i] - lm_pos[ idxs[[i]][j] ] + 1
      
      sum <- 0
      s_num <- 0
      s <- 0
      d <- 0
      I <- 0
      cl <- length(cigar_splited[idxs[[i]][j],])
      for(k in seq(from=2, to=cl, by=2)) {
        if (cigar_splited[idxs[[i]][j] , k] == "") break
        op <- cigar_splited[idxs[[i]][j] , k]
        val <- as.integer(cigar_splited[idxs[[i]][j] , k-1])
        
        if (op == "S") {
          s_num <- s_num + 1
          if (s_num == 2) break
          s <- val
          
        }else if (op == "M") {
          sum <- sum + val
          
        }else if (op == "D") {
          if (pos[j] <= sum) {
            break
          }else {
            d <- d + val
            sum <- sum + val
          }
          
        }else if (op == "I") {
          if (pos[j] <= sum) {
            break
          }else {
            I <- I + val
            sum <- sum + val
          }
          
        }
        
      }
      
      pos[j] <- pos[j] + s - d + I
      qual[j] <- utf8ToInt(substring(as.vector(sam_data[[1]]$qual[j]), pos[j],pos[j]))
    }
    p <- c(p, pos)
    q <- c(q, qual)
  }
}

rm(cigar_splited, lm_pos, pos, qual, rm_pos, vcf, param, reported_pos)
gc(verbose = FALSE, full = TRUE)
t <- as.matrix(sam_data[[1]]$qual[unique(unlist(idxs))])
rm(sam_data)
t <- as.data.frame(t[sample(nrow(t), 100000),])
colnames(t) <- c(1:ncol(t))

t <- melt(t)
colnames(t) <- c("Positions", "Qual")
ggplot(data = t, aes(x=Positions, y=Qual)) +
  geom_boxplot(aes(fill=Positions))+
  xlab("Position")+
  ylab("Qual")+
  labs(title = "Used reads")




# for (i in 1 : length(idxs)) {
#   l <- c(l, lm_pos[idxs[[i]]] - (lm_pos[idxs[[i]]] - 1))
#   pp <- c(pp, reported_pos[i] - (lm_pos[idxs[[i]]] - 1))
#   r <- c(r, rm_pos[idxs[[i]]] - (lm_pos[idxs[[i]]] - 1))
# }

hist(p, plot = T, main = "Histogram of Positions of Q = 30")

hist(q - 33,axis(side = 1, at = seq(1, 75, by = 1), labels = T, tcl = -0.2),  plot = T, main = "Histogram of Qualities of Q = 30")

# hist(p, plot = T)
