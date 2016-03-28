#! /usr/local/bin/Rscript

# ./genepix_two_channel_analysis.R -i ~/workspace/walid/S038/INPUT_FILES -t ~/workspace/walid/S038/S038-targets.txt -l ~/workspace/walid/S038/GQSR02-S1-0001-genepix.gal -o ~/workspace/walid/S038

# http://cran.stat.sfu.ca/
# set a CRAN mirror
local({r <- getOption("repos")
    r["CRAN"] <- "http://cran.stat.sfu.ca/"
    options(repos=r)})

if(!require("plot3D")){
    install.packages("plot3D", dependencies=TRUE);
}

if(!require("getopt")){
    install.packages("getopt", dependencies=TRUE);
}

library("plot3D");
library("getopt");


#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
'input_file','i', 1, "character",
'output_file','o', 1, "character",
'help','h', 0, "logical"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

# if help was asked for print a friendly message
# and exit with a non-zero error code
if(!is.null(opt$help)){
    cat(getopt(spec, usage=TRUE));
    q(status=1);
}

# if no parameter was given for any of these print
# a friendly message and exit with a non-zero error
# code
if(is.null(opt$input_file)){
    cat(getopt(spec, usage=FALSE));
    q(status=1);
}

if(is.null(opt$output_file)){
    cat(getopt(spec, usage=FALSE));
    q(status=1);
}

# Initialize directory and file name variables
input_file <- opt$input_file;
output_file <- opt$output_file;

data = read.table(input_file, header=TRUE, sep="\t", col.names=c("run", "pressure", "temp", "time", "yield", "conc"));

run <- data$run;
pressure <- data$pressure;
temp <- data$temp;
time <- data$time;
yield <- data$yield;
conc <- data$conc;

png(file = output_file, width=800, height=600, res=NA, units="px");
scatter3D(x = temp, y = time, z = pressure, colvar=yield, phi = 0, pch = 20, cex = 2, xlab = "Temperature (°C)", ylab = "Time (min)", zlab="Pressure (bar)", clab = "Yield (g)", main=expression(paste("Supercritical ", CO[2], " Cannabinoid Yield")), labels=run, ticktype = "detailed");
dev.off();
