library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)
q_file = args[[1]]
output_file = args[[2]]
populations = args[[3]]

tbl <- read.table(q_file)

if ( populations >= 3 && populations <= 12 ) {
    colors = brewer.pal(populations, 'Paired')
} else {
    colors = rainbow(populations)
}

pdf(file=output_file, onefile=TRUE, width=7, height=3)
barplot(t(as.matrix(tbl)), col=colors, xlab="Individual #", ylab="Ancestry", border=NA)

dev.off()
