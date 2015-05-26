args <- commandArgs(TRUE);
output_file <- args[1];

x <- read.table('coverage2.txt', skip=1, sep='\t');

individuals <- dim(x)[1];
max_cov <- dim(x)[2] - 2;
max_val <- max(x[-1]) / 100;
colors <- rainbow(individuals);

line_width = 3;
xt = t(x);

xvals <- c(0:max_cov);
values <- as.numeric(as.vector(xt[,1][-1]))/100;

pdf(file=output_file, onefile=TRUE, width=10, height=6);

plot(xvals, values, type='l', ylim=c(0, max_val), xlim=c(0, max_cov), col=colors[1], lwd=line_width, xlab="Coverage", ylab="Proportion");

if (individuals > 1) {
    for (i in 2:individuals) {
        values <- as.numeric(as.vector(xt[,i][-1]))/100;
        lines(xvals, values, col=colors[i], lwd=line_width);
    }
}


names <- as.vector(t(x[1]));
legend(x='topright', legend=names, fill=colors, bty='n');

dev.off();
