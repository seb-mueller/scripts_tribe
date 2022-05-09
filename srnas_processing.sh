cd /projects/TRIBE/srnas/dcl2_zhengming
cd /projects/TRIBE/srnas/trimmed_set3_4/
## counting individual smallRNAs (mapping indepedent). After trimming before mapping:

for file in *_trimmed.fq.gz; do
    zcat $file | grep "^@" -A 1 --no-group-separator | grep -v "^@" | sort | uniq -c > ${file%.fq.gz}.countlist
done

touch tmp3
for file in 01*countlist; do
  echo $file
  cp tmp3 tmp4
  join -j2 tmp5 $file -a 1 -a 2 -e '0' -o '0,1.1,2.1' > tmp3
  sort tmp4 > tmp5
done

  join -j2 $file tmp3 -a 1 -a 2 -e '0' -o '0,1.1,2.1' > tmp3

https://www.theunixschool.com/2012/01/join-command.html?m=1

##in R
library(stringr)
library(dplyr)
library(purrr)
library(Biostrings)
path="/projects/TRIBE/srnas/dcl2_zhengming/"
path="/projects/TRIBE/srnas/trimmed_set3_4/"
fls.countlist <- list.files(path, pattern="countlist$", full.names =F)

idx=which(str_detect(fls.countlist, ".*KP[:digit:]{4}.*"))  # finds +DESL -DESL
fls.countlist =fls.countlist[c(1:4, idx)]

# half
idx=idx[c(TRUE,FALSE)]
fls.countlist=fls.countlist[c(1,3, idx)]

map(fls_list_counts, head)
tmp2 <- Reduce(merging2,tmp)     #merge into a count_table for all sRNA species

merging<-function(x,y) merge(x,y,by='sequence',all=T)
merging<-function(x,y) merge(x,y,by='V2',all=T)
merging2<-function(x,y) left_join(x,y,by='sequence',all=T)
fComb <- function(y,x) left_join(y,fRead(x))
df    <-  Reduce(fComb,myfiles[-1],fRead(myfiles[1]))

fls_list_counts <- list()
for (file in fls.countlist) { 
  # file=fls.countlist[1]
  message(file)
  tmp <- read.table(file.path(path, file), row.names=NULL,colClasses = c("numeric","character")) 
  tmp3=width(tmp$V2)
  tmp2 <- tmp[tmp3>19 & tmp3 < 26,]
  # table(width(tmp2$V2))
  colnames(tmp2) <- c(file,  "sequence")
  fls_list_counts[[file]] <- tmp2
}

map(fls_list_counts, table)

cnt <- Reduce(merging2,fls_list_counts)     #merge into a count_table for all sRNA species
cnt <- Reduce(merging,fls_list_counts)     #merge into a count_table for all sRNA species
colnames(cnt) <- c("sequence", names(fls_list_counts))
colnames(cnt) <- str_replace(colnames(cnt), "_pool.*list", "")
colnames(cnt) <- str_replace(colnames(cnt), "_trim.*list", "")
tmp2 <- merging(fls_list_counts[[1]],fls_list_counts[[2]])     #merge into a count_table for all sRNA species
map(fls_list_counts,head)     #merge into a count_table for all sRNA species
cnt <- Reduce(merging,tmp)     #merge into a count_table for all sRNA species
library(data.table)
tmp1=data.table(tmp[[1]])
tmp2=data.table(tmp[[2]])
tmp1[tmp2, on = 'sequence', bb := i.]

cnt$width <- width(cnt$sequence)

cnt_dcl2 <- cnt %>%
   replace(is.na(.), 0) %>%
   relocate(sequence, width)

save(cnt,file="cnt_tribe.rdata")
save(cnt_dcl2,file="cnt_dcl2.rdata")
# [1] "/projects/TRIBE/srnas/dcl2_zhengming"


tmp <- colnames(cnt)
colnames(cnt)=str_replace(colnames(cnt), "-.*","")
colnames(cnt)=str_replace(colnames(cnt), ".*_","")
#colnames(cnt) <- c("sequence",paste(sapply(strsplit(fls.countlist,"_"),function(x) paste(x[1],x[2],x[3],sep="_")),sep="_"))

tmp=cnt %>%
  group_by(width) %>%
  summarize(across(2:23, sum))
  # summarize(sum = sum(KP4041))
colnames(tmp)=str_replace(colnames(tmp), "[0-9]*_","")
colnames(tmp)=str_replace(colnames(tmp), "_.*","")
# colnames(tmp)=colnames(cnt)

cnt_summarized <- tmp %>%
      pivot_longer(
                   cols = starts_with("K"),
                   names_to = c("plant","Rep"),
                   names_pattern = "K(.*)-(.*)",
                   values_to = "count",
                   values_drop_na = TRUE,
      )
cnt_summarized$size <- as.factor(cnt_summarized$sequence)

colnames(cnt_summarized) <- c("sequence",groupdf.unique$Library)

save(cnt_summarized,file="cnt_summarized.rdata")
  # sequence plant  Rep     count size
# 1       20 M82C   1      659384 20
# 2       20 M82C   2      786459 20
# 3       20 PenneC 1      637783 20
# 5       20 P1511  1     1199087 20


colnames(cnt) <- c("sequence",groupdf.unique$Library)
cnt[is.na(cnt)]<-0
cnt$width <- nchar(cnt$sequence)
rownames(cnt) <- cnt$sequence
cnt <- cnt[,-1]
cnt.high <- cnt[rowSums(cnt[,1:24])>100,] 
##annoating miRNAs
#as.character(microRNAs) %in% as.character(microRNAs)[1]
cnt$is.miRNA <- rownames(cnt) %in% as.character(microRNAs)
cnt$miRNA <- cnt$is.miRNA 
for (mi in which(cnt$is.miRNA)) {
	cnt$miRNA[mi] <- names(microRNAs)[as.character(microRNAs) %in% rownames(cnt)[mi]][1]
}
##looking at data
cnt_agg <- aggregate(cnt[,1:24], list(size = cnt$width), sum)
cnt_agg_pcts <-  t(t(cnt_agg[,-1]) * 100 / colSums(cnt_agg[,-1]))
cnt_agg_pcts2025 <-  t(t(cnt_agg[cnt_agg$size %in% 20:25,-1]) * 100 / colSums(cnt_agg[cnt_agg$size %in% 20:25,-1]))
rownames(cnt_agg_pcts2025) <- 20:25


