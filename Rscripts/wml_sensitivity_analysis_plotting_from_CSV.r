#!/usr/bin/env Rscript

#
# This script plots a boxchart for the provided CSV data.
# It considers each input file as one group that will be contrasted to
# the others in the scatterplot.
# All input CSV-files must have the same header.
#

library(ggplot2)
library(reshape2)
library(data.table)
library(methods)
library(car)
library(coin)
library(stats)
library(sjstats) # anova effect size "anova_stats" function
library(rstatix) # kruskall wallis effect size
library(effsize)
library(rcompanion)

#
# This class reads in a sequence of CSV files, stacks the included data
# and allows to trim the columns according to a provided list of desired
# columns to be included in the trimmed output. The data from each CSV
# file is treated as one group, i.e. data wihtin a CSV file comes from
# subjects of the same group, mulitple CSV files represent multipe groups
# of subjects.
#
# Instantiation requires the list of (paths to) the CSV files and a flag
# denoting whether the data should be treated as groups (one CSV = one group).
#
DataProvider <- setRefClass(
    "DataProvider",

    fields = list(
        input_CSVs      = "character",
        stacked_input   = "data.frame",
        trimmed_input   = "data.frame",
        data_is_grouped = "logical"
    ),

    methods = list(
        # initialize method acts as the constructor
        initialize = function( list_of_input_CSVs, is_grouped ) {
            input_CSVs <<- list_of_input_CSVs;

            input <- list()
            for( i in seq(input_CSVs) ){
                input[[i]] <- tryCatch(read.csv( file=input_CSVs[[i]] ), error=function(e){print("File does not exist! Skipping input file")})

                # add a group column to this dataframe to identify which
                # group it is coming from
                if( is_grouped ) {
                    input[[i]]$group <- as.factor(i-1)
                }
            }

            stacked_input   <<- rbindlist( input, use.names=TRUE, fill=TRUE );
            trimmed_input   <<- stacked_input;
            data_is_grouped <<- is_grouped;
        },

        # This function removes all unwanted columns from the input data
        # frame and saves the result in 'trimmed_input'
        trimStackedColumns = function( cols_to_be_used ) {
            trimmed_input <<- stacked_input;

            if( ! missing(cols_to_be_used) ) {
                if( data_is_grouped ){  # always keep the group if present
                    cols_to_be_used = c( cols_to_be_used, "group" )
                }
                for( col in setdiff( colnames(trimmed_input), c( cols_to_be_used ) ) ){
                    trimmed_input[ which(names(trimmed_input) == col) ] <<- NULL;
                }
            }
        }
    )
)

#
# This custom functions prints equidistant ticks at log-scales.
#
base_breaks <- function(n = 10){
    function(x) {   # 'x' will be fed in by the 'scale'-function and is a 2-element list containing the 'min' and 'max'
        min_max = range(x, na.rm=TRUE)
        min_max = log10(min_max)
        axisTicks(min_max, log = TRUE, n = n)
    }
}

#
# This functions plots a boxplot with the individual data scattered
# along the box.
#
plot = function(
    data,
    use_log_scale,
    main_title,
    sub_title,
    legend_title,
    y_axis_title,
    x_axis_title,
    output_path,
    file_name )
{
    scale_function <- scale_y_continuous

    if( use_log_scale ){
        #scale_function <- function(){ scale_y_log10(name = y_axis_title,labels = scales::comma, breaks=base_breaks()) }
        scale_function <- function(){ scale_y_log10(name = y_axis_title,breaks=base_breaks()) }
    }

    # use the following 4 lines to create a common lower boundary (1e-9) for the "all-tissues"-Sobl inxed plots
    group<-rep(0,5)
    variable<-as.factor(c('Skin', 'Skull', 'Cerebrospinal fluid', 'Gray matter', 'White matter'))
    value<-rep(1.8e-9, 5)
    dummy_data<-data.frame(group, variable, value)

    ggplot( data, aes( x=variable, y=value, fill=factor(group), ordered=TRUE) ) +
    geom_boxplot(alpha=.9, outlier.shape = 32) +
    scale_function()+
    labs(fill = "group", title=main_title, subtitle=sub_title )+
    guides(fill=guide_legend(title=legend_title))+
    geom_point(position=position_jitterdodge(jitter.width=.25), alpha=0.3, data=data[ with(data, ! is.na(value)), ] ) +
    geom_point(data=dummy_data, aes(x=variable, y=value), color="darkorange3", alpha=0, guide=FALSE, position=position_jitterdodge(jitter.width=.25)) +
    theme_bw() +
    theme(plot.title = element_text(lineheight=.8, size = 18, face="bold", hjust = 0.5),axis.title=element_text(size=14, face="italic"))+
    annotation_logticks(sides="l") +
    scale_x_discrete(name = x_axis_title) +
    expand_limits( y = c(0,0.002) ) +
    scale_fill_brewer(palette = "Greens")

    ggsave(file_name, plot=last_plot(), device="pdf", path=output_path, dpi=300)
}

#
# This functions renames the columns "lesion", "csf", "skull", "scalp",
# "wm", "gm", "var", "mean" of a data frame to more descriptive names
#
rename_columns = function( data_frame ) {
    colnames(data_frame)[which(names(data_frame) == "lesion")] <- "Lesion"
    colnames(data_frame)[which(names(data_frame) == "csf")] <- "Cerebrospinal fluid"
    colnames(data_frame)[which(names(data_frame) == "skull")] <- "Skull"
    colnames(data_frame)[which(names(data_frame) == "scalp")] <- "Skin"
    colnames(data_frame)[which(names(data_frame) == "wm")] <- "White matter"
    colnames(data_frame)[which(names(data_frame) == "gm")] <- "Gray matter"
    colnames(data_frame)[which(names(data_frame) == "var")] <- "Variance"
    colnames(data_frame)[which(names(data_frame) == "mean")] <- "Mean field strength"

    return(data_frame)
}

########
# MAIN #
########
DISPLAY_IN_PERCENT=FALSE
OUTPUT_DIR="/temp/"
USE_LOG_SCALE=TRUE

args = commandArgs(trailingOnly=TRUE)

if( length(args) == 0 )
{
    stop( "At least one CSV file with data to be plotted must be provided.", call.=FALSE );
}

reader <- DataProvider( args, TRUE );
all_data <- reader$stacked_input
all_data$subjectID <- NULL
if(FALSE){
## show means and variance and standard deviations
cat( "##### DESCRIPTIVE STATISTICS #####\n" )
for(fazekas_score in c(0,1,2,3)){
    cat("*** Fazekas ",fazekas_score," ***\n")
    cat(" - Electrical field strength\n")
    cat("\tmean=",mean(all_data$mean[all_data$group == fazekas_score]),"\n")
    cat("\tvariance=",var(all_data$mean[all_data$group == fazekas_score]),"\n")
    cat("\tstandard deviation=",sd(all_data$mean[all_data$group == fazekas_score]),"\n")
    cat(" - Variance of electrical field strength\n")
    cat("\tmean=",mean(all_data$var[all_data$group == fazekas_score]),"\n")
    cat("\tvariance=",var(all_data$var[all_data$group == fazekas_score]),"\n")
    cat("\tstandard deviation=",sd(all_data$var[all_data$group == fazekas_score]),"\n")
    cat(" - Normalized sobol index scalp\n")
    cat("\tmean=",mean(all_data$scalp[all_data$group == fazekas_score]),"\n")
    cat("\tvariance=",var(all_data$scalp[all_data$group == fazekas_score]),"\n")
    cat("\tstandard deviation=",sd(all_data$scalp[all_data$group == fazekas_score]),"\n")
    cat(" - Normalized sobol index skull\n")
    cat("\tmean=",mean(all_data$skull[all_data$group == fazekas_score]),"\n")
    cat("\tvariance=",var(all_data$skull[all_data$group == fazekas_score]),"\n")
    cat("\tstandard deviation=",sd(all_data$skull[all_data$group == fazekas_score]),"\n")
    cat(" - Normalized sobol index csf\n")
    cat("\tmean=",mean(all_data$csf[all_data$group == fazekas_score]),"\n")
    cat("\tvariance=",var(all_data$csf[all_data$group == fazekas_score]),"\n")
    cat("\tstandard deviation=",sd(all_data$csf[all_data$group == fazekas_score]),"\n")
    cat(" - Normalized sobol index gm\n")
    cat("\tmean=",mean(all_data$gm[all_data$group == fazekas_score]),"\n")
    cat("\tvariance=",var(all_data$gm[all_data$group == fazekas_score]),"\n")
    cat("\tstandard deviation=",sd(all_data$gm[all_data$group == fazekas_score]),"\n")
    cat(" - Normalized sobol index wm\n")
    cat("\tmean=",mean(all_data$wm[all_data$group == fazekas_score]),"\n")
    cat("\tvariance=",var(all_data$wm[all_data$group == fazekas_score]),"\n")
    cat("\tstandard deviation=",sd(all_data$wm[all_data$group == fazekas_score]),"\n")
    cat(" - Normalized sobol index lesion\n")
    cat("\tmean=",mean(all_data$lesion[all_data$group == fazekas_score]),"\n")
    cat("\tvariance=",var(all_data$lesion[all_data$group == fazekas_score]),"\n")
    cat("\tstandard deviation=",sd(all_data$lesion[all_data$group == fazekas_score]),"\n")
}

## tests ##
cat( "##### SIGNIFICANCE TESTS #####" )
# normality
sapply( all_data, function( col ) shapiro.test( as.numeric(col)) )
for(col in setdiff(names(all_data), "group")){
    pdf(paste( OUTPUT_DIR,"qqplot_",col,".pdf", sep=""));
    qqnorm(all_data[[col]]);
    dev.off();
}

# homogeneity of variance
apply(all_data, 2, function(x) leveneTest( y = as.numeric(x), group = all_data$group ) )

# parametric test: ANOVA
lapply(setdiff(names(all_data), "group"), function( col ) {
    print(paste0( "Dataset: ", col ));
    fit<-aov( as.formula( paste( col, " ~ group") ), all_data);
    print(anova_stats(fit)) # explicitly printing this causes a double console output - no big deal
    print(cohens_d(all_data, as.formula( paste( col, " ~ group")), paired=TRUE ))
})

# non-parametric test: Kruskal–Wallis one-way analysis of variance
# as an alternative to the Mann-Whitney-U test with more than two groups
lapply(setdiff(names(all_data), "group"), function( col ) {
    print(paste0( "Dataset: ", col ));
    print(kruskal.test( as.formula( paste( col, " ~ group") ), all_data))
    print(kruskal_effsize(all_data,as.formula( paste( col, " ~ group")))) # explicitly printing this causes a double console output - no big deal
})

############### Early quitting to skip plotting!!! #############################
#quit()
}
## plotting ##
# [1] contributions of all tissues
reader$trimStackedColumns( c("scalp","skull","csf","gm","wm","lesion") );
to_be_plotted <- reader$trimmed_input;

if( DISPLAY_IN_PERCENT ) {
    to_be_plotted$scalp<-to_be_plotted$scalp*100
    to_be_plotted$skull<-to_be_plotted$skull*100
    to_be_plotted$csf<-to_be_plotted$csf*100
    to_be_plotted$gm<-to_be_plotted$gm*100
    to_be_plotted$wm<-to_be_plotted$wm*100
    to_be_plotted$lesion<-to_be_plotted$lesion*100
}

to_be_plotted <- rename_columns( to_be_plotted );
# melt the dataframe down to 3 columns:
#   1) group-colum,
#   2) variable-name column,
#   3) data column
# => Data from multiple variables will be stacked.
to_be_plotted=melt(to_be_plotted,id.vars='group');

ordering_row_idx<-c(    # create a list of indices in the order that we want
    c(which(to_be_plotted$variable == "Skin")),
    c(which(to_be_plotted$variable == "Skull")),
    c(which(to_be_plotted$variable == "Cerebrospinal fluid")),
    c(which(to_be_plotted$variable == "Gray matter")),
    c(which(to_be_plotted$variable == "White matter")),
    c(which(to_be_plotted$variable == "Lesion"))
)
to_be_plotted<-to_be_plotted[ordering_row_idx,]
#to_be_plotted<-to_be_plotted[order(as.character(to_be_plotted$variable)),] # for alphabetical ordering
to_be_plotted$variable <- factor(to_be_plotted$variable, levels=unique(to_be_plotted$variable)) # must be marked as factor so ggplot will order it alphabetically...


plot(
    data = to_be_plotted,
    use_log_scale = USE_LOG_SCALE,
    main_title = "Contribution of individual tissue to total variance",
    sub_title = "Groups: Fazekas 0 vs. 1 vs. 2 vs. 3",
    legend_title = "Fazekas\nscore",
    y_axis_title = "Share of total variance",
    x_axis_title = "Sobol index",
    output_path = OUTPUT_DIR,
    file_name = "contributions_of_all_tissues.pdf"
);

# [2] contributions of scalp only
reader$trimStackedColumns( c("scalp") );
to_be_plotted <- reader$trimmed_input;
to_be_plotted <- rename_columns( to_be_plotted );
to_be_plotted=melt(to_be_plotted,id.vars='group');

plot(
    data = to_be_plotted,
    use_log_scale = USE_LOG_SCALE,
    main_title = "Contribution of scalp tissue to total variance",
    sub_title = "Groups: Fazekas 0 vs. 1 vs. 2 vs. 3",
    legend_title = "Fazekas\nscore",
    y_axis_title = "Share of total variance",
    x_axis_title = "Sobol index",
    output_path = OUTPUT_DIR,
    file_name = "contributions_of_scalp_tissues.pdf"
);

# [3] contributions of skull only
reader$trimStackedColumns( c("skull") );
to_be_plotted <- reader$trimmed_input;
to_be_plotted <- rename_columns( to_be_plotted );
to_be_plotted=melt(to_be_plotted,id.vars='group');

plot(
    data = to_be_plotted,
    use_log_scale = USE_LOG_SCALE,
    main_title = "Contribution of skull tissue to total variance",
    sub_title = "Groups: Fazekas 1 vs. 2 vs. 3",
    legend_title = "Fazekas\nscore",
    y_axis_title = "Share of total variance",
    x_axis_title = "Sobol index",
    output_path = OUTPUT_DIR,
    file_name = "contributions_of_skull_tissues.pdf"
);

# [4] contributions of CSF only
reader$trimStackedColumns( c("csf") );
to_be_plotted <- reader$trimmed_input;
to_be_plotted <- rename_columns( to_be_plotted );
to_be_plotted=melt(to_be_plotted,id.vars='group');

plot(
    data = to_be_plotted,
    use_log_scale = USE_LOG_SCALE,
    main_title = "Contribution of cerebrospinal-fluid to total variance",
    sub_title = "Groups: Fazekas 0 vs. 1 vs. 2 vs. 3",
    legend_title = "Fazekas\nscore",
    y_axis_title = "Share of total variance",
    x_axis_title = "Sobol index",
    output_path = OUTPUT_DIR,
    file_name = "contributions_of_csf_tissues.pdf"
);

# [5] contributions of GM only
reader$trimStackedColumns( c("gm") );
to_be_plotted <- reader$trimmed_input;
to_be_plotted <- rename_columns( to_be_plotted );
to_be_plotted=melt(to_be_plotted,id.vars='group');

plot(
    data = to_be_plotted,
    use_log_scale = USE_LOG_SCALE,
    main_title = "Contribution of gray matter tissue to total variance",
    sub_title = "Groups: Fazekas 0 vs. 1 vs. 2 vs. 3",
    legend_title = "Fazekas\nscore",
    y_axis_title = "Share of total variance",
    x_axis_title = "Sobol index",
    output_path = OUTPUT_DIR,
    file_name = "contributions_of_gm_tissues.pdf"
);

# [6] contributions of WM
reader$trimStackedColumns( c("wm") );
to_be_plotted <- reader$trimmed_input;
to_be_plotted <- rename_columns( to_be_plotted );
to_be_plotted=melt(to_be_plotted,id.vars='group');

plot(
    data = to_be_plotted,
    use_log_scale = USE_LOG_SCALE,
    main_title = "Contribution of whate matter tissue to total variance",
    sub_title = "Groups: Fazekas 0 vs. 1 vs. 2 vs. 3",
    legend_title = "Fazekas\nscore",
    y_axis_title = "Share of total variance",
    x_axis_title = "Sobol index",
    output_path = OUTPUT_DIR,
    file_name = "contributions_of_wm_tissues.pdf"
);

# [7] contributions of WMLs only
reader$trimStackedColumns( c("lesion") );
to_be_plotted <- reader$trimmed_input;
to_be_plotted <- rename_columns( to_be_plotted );
to_be_plotted=melt(to_be_plotted,id.vars='group');

plot(
    data = to_be_plotted,
    use_log_scale = USE_LOG_SCALE,
    main_title = "Contribution of lesioned tissue to total variance",
    sub_title = "Groups: 1 vs. 2 vs. 3",
    legend_title = "Fazekas\nscore",
    y_axis_title = "Share of total variance",
    x_axis_title = "Sobol index",
    output_path = OUTPUT_DIR,
    file_name = "contributions_of_lesioned_tissues.pdf"
);

# [8] mean electrical field strength
reader$trimStackedColumns( c("mean") );
to_be_plotted <- reader$trimmed_input;
to_be_plotted <- rename_columns( to_be_plotted );
to_be_plotted=melt(to_be_plotted,id.vars='group');

plot(
    data = to_be_plotted,
    use_log_scale = USE_LOG_SCALE,
    main_title = "Mean electrical field strength in ROI",
    sub_title = "Groups: Fazekas 0 vs. 1 vs. 2 vs. 3",
    legend_title = "Fazekas\nscore",
    y_axis_title = "Electrical field strength",
    x_axis_title = "",
    output_path = OUTPUT_DIR,
    file_name = "electrical_field_strength.pdf"
);

# [9] variance of electrical field strength
reader$trimStackedColumns( c("var") );
to_be_plotted <- reader$trimmed_input;
to_be_plotted <- rename_columns( to_be_plotted );
to_be_plotted=melt(to_be_plotted,id.vars='group');

plot(
    data = to_be_plotted,
    use_log_scale = USE_LOG_SCALE,
    main_title = "Variance of the electrical field strength in ROI",
    sub_title = "Groups: Fazekas 0 vs. 1 vs. 2 vs. 3",
    legend_title = "Fazekas\nscore",
    y_axis_title = "Variance of the electrical field strength",
    x_axis_title = "",
    output_path = OUTPUT_DIR,
    file_name = "variance.pdf"
);
