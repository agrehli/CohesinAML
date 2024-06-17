#!/bin/bash
# by Michael Rehli, Mai 13, 2018
# bash script to generate a  hist plot (pdf) including 95% confidence band
# fixed issue with neg. coordinates, Nov. 2021
# 

#setting PATH
PATH_R=${DIR_SOFT}/package/RBioC/stretch_R-4.0.3_Bioc-3.12_intel-mkl-64bit-2020.1-102/lib/R/bin
export PATH=${PATH_R}:${PATH}
export PATH

TMPDIR="/loctmp"

# required input:
# 
# -g <homer derived ghist file>
# -w <plot width> (pixels for png)
# -h <plot height> (pixels for png)
# -c <color> (default "blue")
# -m <max value for color scale> (default 10)
# -d <output directory>
# -n <name>

# Set Script Name variable
SCRIPT=`basename ${BASH_SOURCE[0]}`
calc() { awk "BEGIN{print $*}"; }

# Set fonts for Help.
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`

# USAGE function
function USAGE {
  echo -e \\n"Help documentation for ${BOLD}${SCRIPT}.${NORM}"\\n
  echo -e "${REV}Basic usage:${NORM}\n${BOLD}$SCRIPT -g <homer derived ghist file(s)> -s <sample name(s)> -w <individual plot width> -d <output directory> -n <name> -c <color(s)> -x <window> -y <y axis range> ${NORM}"\\n
  echo "Required:"
  echo "${BOLD}-g${NORM} <homer derived ghist file or list of files (space separated and in quotation marks)>"
  echo "${BOLD}-s${NORM} <sample name or list of sample names (space separated and framed by quotation marks)>"
  echo "${BOLD}-d${NORM} <output directory>"
  echo "${BOLD}-n${NORM} <name> used for output files"
  echo "${BOLD}-c${NORM} <color or list of colors (space separated, one for each sample, framed by quotation marks)> (default blue)"
  echo "${BOLD}-x${NORM} <X axis window size (in bp)>"
  echo "Optional:"
  echo "${BOLD}-m${NORM} <move center (in +/- bp)>"
  echo "${BOLD}-w${NORM} <plot width (in inch, default: 2)>"
  echo "${BOLD}-h${NORM} <plot hight (in inch, default: 2)>"
  echo "${BOLD}-b${NORM} <line width (default: .5)>"
  echo "${BOLD}-t${NORM} <legend text size (default: 2)>"
  echo "${BOLD}-u${NORM} <method> for fit curve (default: loess), optional: glm, lm, gam, rlm >"
  echo "${BOLD}-l${NORM} <span> for fit curve (default: 0.04)>"
  echo "${BOLD}-f${NORM} only print fit curve>"
  echo "${BOLD}-a${NORM} <x pos.> add vertical line(s) at position (default: none)>"
  echo "${BOLD}-e${NORM} <color(s)> add vertical line(s) at position (default: none)>"
  echo -e "${BOLD}-y${NORM} <Y axis range (min and max values used for Y axis, space separated and framed by quotation marks)>"\\n
  echo -e "Example: ${BOLD}$SCRIPT -f testing/ghist_ATAC.txt -d /misc/data/analysis/testing -n test1 -c firebrick -x 1000 -y \"0 15\" ${NORM}"\\n
  exit 1
}

#Check the number of arguments. If none are passed, print help and exit.
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
  USAGE
fi

# Set defaultS
COLOR="BLUE"
MAXVALUE=10
GIVEN=0
SHIFT=0
PLOTHIGHT=2
PLOTWIDTH=2
LINEWIDTH=.5
LABELSIZE=2
FIT=0
FITONLY=0
SPAN=0.04
VLINE=0
METHOD="loess"

# Parse command line options
while getopts :g:s:d:n:c:x:m:h:w:y:b:t:fl:a:e: OPTIONS; do
  case $OPTIONS in
    g)  #set option "f" 
      GHIST=($OPTARG)
      SAMPLENUMBER=${#GHIST[@]}
      echo ${#GHIST[@]}
      ;;
    s)  #set option "s"
      SAMPLES=($OPTARG)
      ;;
    d)  #set option "d"
      DIRECTORY=$OPTARG
      ;;
    n)  #set option "n"
      NAME=$OPTARG
      ;;
    c)  #set option "c"
      COLORS=($OPTARG)
      ;;
    y)  #set option "y"
      AXISVALUES=($OPTARG)
#	  MIDVALUED=$MAXVALUE/2
#	  MIDVALUEU=$MAXVALUE/2+0.01
	;;
    x)  #set option "x"
      WINDOW=$OPTARG
    ;;
    m)  #set option "m"
      SHIFT=$OPTARG
	;;
    h)  #set option "h"
      PLOTHIGHT=$OPTARG
	;;
    w)  #set option "w"
      PLOTWIDTH=$OPTARG
	;;
    t)  #set option "t"
      LABELSIZE=$OPTARG
	;;
    b)  #set option "b"
      LINEWIDTH=$OPTARG
	;;
    l)  #set option "l"
      SPAN=$OPTARG
      FIT=1
	;;
    f)  #set option "f"
      FITONLY=1
	;;
    a)  #set option "a"
      VLINEPOS=($OPTARG)
      VLINENUMBER=${#VLINEPOS[@]}
      VLINE=1
	;;
    e)  #set option "e"
      VLINECOLORS=($OPTARG)
	;;
    \?) #unrecognized option - show USAGE
      echo -e \\n"Option -${BOLD}$OPTARG${NORM} not allowed."
      USAGE
      ;;
  esac
done
shift $((OPTIND-1)) 

# Sanity checks and definition of variables needed for R sript

# Check number of samples

if [ ${#GHIST[@]} != ${#COLORS[@]} ] ; then
	echo "The number of colors given (${#COLOR[@]}) doesn't match the sample number! (${#GHIST[@]})"
    echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
    exit 3
fi

# if [ ${#GHIST[@]} != ${#MAXVALUE[@]} ] ; then
# 	echo "The number of max. colors values given (${#MAXVALUE[@]) doesn't match the sample number! (${#GHIST[@])"
#     echo -e "Use ${BOLD}$SCRIPT -h${NORM} to see the help documentation."\\n
#     exit 3
# fi

# Define axis limits

YAXISD=${AXISVALUES[0]}
YAXISU=${AXISVALUES[1]}

XAXISMAXD=$(calc "(${WINDOW}/-2)-($SHIFT)")
XAXISMAXU=$(calc "(${WINDOW}/2)-($SHIFT)")
ANNSTEP=$(calc "((${AXISVALUES[1]})-(${AXISVALUES[0]}))/11")
echo $ANNSTEP

echo "Generating R Script"

# Generate and execute R script

cat >"${DIRECTORY}/R.hist.${NAME}.R" <<EOF
setwd("${DIRECTORY}")
library(ggplot2)
data0 <- read.table("${GHIST[0]}", header=T, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)
tdata0 <- t(data0[,-1])
m0 <- apply(tdata0, 1, mean)
s <- apply(tdata0, 1, sd)
n <- nrow(as.matrix(tdata0))
###### calculating the confidence interval
error <- qnorm(0.975)*s/sqrt(n)
###### upper and lower boundaries
ci975 <- m0+error
ci025 <- m0-error
d0 <- as.numeric(as.character(rownames(tdata0)))
result0 <- data.frame(cbind(d0,m0, s,ci975,ci025))
#head(result0)
#fit <- loess(m0~d0,data=result0,family = "gaussian",span = 0.02)
#result0$pred <- predict(fit)

###### starting the plot
p <- ggplot(result0, aes(d0,m0))
p <- p + theme_bw(base_size=12) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
EOF

if [ ${VLINE} == 1 ] ; then
for (( i = 0 ; i < ${VLINENUMBER} ; i++ )); do
cat >>"${DIRECTORY}/R.hist.${NAME}.R" <<EOF
p <- p + geom_vline(xintercept=${VLINEPOS[$i]}, linetype="dashed", color = "${VLINECOLORS[$i]}", size = .25)
EOF
done
fi

if [ ${FITONLY} == 0 ] ; then
cat >>"${DIRECTORY}/R.hist.${NAME}.R" <<EOF
p <- p + geom_ribbon(aes(x=d0, ymax=ci975, ymin=ci025), fill="${COLORS[0]}", alpha=.15)
EOF
fi

if [ ${FIT} == 0 ] ; then
cat >>"${DIRECTORY}/R.hist.${NAME}.R" <<EOF
p <- p + geom_line(colour="${COLORS[0]}", size=${LINEWIDTH}) 
EOF
else
cat >>"${DIRECTORY}/R.hist.${NAME}.R" <<EOF
p <- p + stat_smooth(method = "${METHOD}", geom = "smooth", span = ${SPAN}, n = ${WINDOW}, se = TRUE, fullrange = TRUE, colour="${COLORS[0]}", fill = "${COLORS[0]}", size = ${LINEWIDTH}) 
EOF
fi
cat >>"${DIRECTORY}/R.hist.${NAME}.R" <<EOF
p <- p + xlab("Distance to center") + ylab("Coverage")
# p <- p + xlim(${XAXISMAXD}, ${XAXISMAXU}) + ylim(${YAXISD}, ${YAXISU}) 
p <- p + scale_x_continuous(expand=c(0.01,0.01), limits = c(${XAXISMAXD}, ${XAXISMAXU})) + ylim(${YAXISD}, ${YAXISU}) 
p <- p + annotate("text", x = ${XAXISMAXD} , y = ${YAXISU}, label ="${SAMPLES[0]}", size = ${LABELSIZE}, colour="${COLORS[0]}", hjust=0, vjust=1)
EOF

if [ ${SAMPLENUMBER} >1 ] ; then
for (( i = 1 ; i < ${SAMPLENUMBER} ; i++ )); do
YANN=$(calc "(${AXISVALUES[1]})-(${ANNSTEP}*${i})")
cat >>"${DIRECTORY}/R.hist.${NAME}.R" <<EOF
data${i} <- read.table("${GHIST[$i]}", header=T, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)
tdata${i} <- t(data${i}[,-1])
m${i} <- apply(tdata${i}, 1, mean)
s <- apply(tdata${i}, 1, sd)
n <- nrow(as.matrix(tdata${i}))
###### calculating the confidence interval
error <- qnorm(0.975)*s/sqrt(n)
###### upper and lower boundaries
ci975 <- m${i}+error
ci025 <- m${i}-error
d${i} <- as.numeric(as.character(rownames(tdata${i})))
result${i} <- data.frame(cbind(d${i},m${i}, s,ci975,ci025))
###### continuing the plot
EOF

if [ ${FITONLY} == 0 ] ; then
cat >>"${DIRECTORY}/R.hist.${NAME}.R" <<EOF
p <- p + geom_ribbon(data=result${i}, aes(x=d${i}, ymax=ci975, ymin=ci025), fill="${COLORS[$i]}", alpha=.15)
EOF
fi

if [ ${FIT} == 0 ] ; then
cat >>"${DIRECTORY}/R.hist.${NAME}.R" <<EOF
p <- p + geom_line(data=result${i}, aes(x=d${i}, y=m${i}), colour="${COLORS[$i]}", size=${LINEWIDTH})
EOF
else
cat >>"${DIRECTORY}/R.hist.${NAME}.R" <<EOF
p <- p + stat_smooth(data=result${i}, aes(x=d${i}, y=m${i}), geom = "smooth", method = "${METHOD}", span = ${SPAN}, n = ${WINDOW}, se = TRUE, fullrange = TRUE, colour="${COLORS[$i]}", fill = "${COLORS[$i]}", size=${LINEWIDTH}) 
EOF
fi

cat >>"${DIRECTORY}/R.hist.${NAME}.R" <<EOF
p <- p + annotate("text", x = ${XAXISMAXD} , y = ${YANN}, label ="${SAMPLES[$i]}", size = ${LABELSIZE}, colour="${COLORS[$i]}", hjust=0, vjust=1)
EOF
done
fi


cat >>"${DIRECTORY}/R.hist.${NAME}.R" <<EOF
pdf(file="${DIRECTORY}/${NAME}.hist.pdf", height=${PLOTHIGHT}, width=${PLOTWIDTH})
plot(p)
dev.off()
EOF
chmod 750 "${DIRECTORY}/R.hist.${NAME}.R"
R < "${DIRECTORY}/R.hist.${NAME}.R"  --no-save 
