library(reshape2)
library(data.table)
library(methods)
library(car)
library(coin)
library(stats)
library(dunn.test)

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
            for( i in seq_along(input_CSVs) ){
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

#~ args=c("/media/simulations/WML/analysis/CSV/InnerROI/FAZEKAS1_inner_ROI_stats_summary.csv","/media/simulations/WML/analysis/CSV/InnerROI/FAZEKAS2_inner_ROI_stats_summary.csv","/media/simulations/WML/analysis/CSV/InnerROI/FAZEKAS3_inner_ROI_stats_summary.csv")
# args=c('/media/simulations/WML/analysis/CSV/wholebrain/FAZEKAS0_wholebrain_stats_summary.csv','/media/simulations/WML/analysis/CSV/wholebrain/FAZEKAS1_wholebrain_stats_summary.csv','/media/simulations/WML/analysis/CSV/wholebrain/FAZEKAS2_wholebrain_stats_summary.csv','/media/simulations/WML/analysis/CSV/wholebrain/FAZEKAS3_wholebrain_stats_summary.csv')
args=c('/media/simulations/WML/analysis/CSV/non-normalized_sobol_indices/FAZEKAS0_wholebrain_stats_summary.csv','/media/simulations/WML/analysis/CSV/non-normalized_sobol_indices/FAZEKAS1_wholebrain_stats_summary.csv','/media/simulations/WML/analysis/CSV/non-normalized_sobol_indices/FAZEKAS2_wholebrain_stats_summary.csv','/media/simulations/WML/analysis/CSV/non-normalized_sobol_indices/FAZEKAS3_wholebrain_stats_summary.csv')
reader<-DataProvider(args,TRUE)
all_data<-reader$stacked_input
all_data$subjectID<-NULL
#~ pairwise.t.test(all_data$skull,all_data$group, p.adj="bonf")     # normalized sobol indices
#~ pairwise.t.test(all_data$gm,all_data$group, p.adj="bonf")        # normalized sobol indices
print("wholebrain: wm")
dunn.test(all_data$wm,all_data$group, method="bonferroni")   # normalized sobol indices
print("wholebrain: lesion")
dunn.test(all_data$lesion,all_data$group, method="bonferroni")   # normalized sobol indices

#~ print("wholebrain: wm")
#~ dunn.test(all_data$wm,all_data$group, method="bonferroni")      # non-normalized sobol indices
#~ print("wholebrain: lesion")
#~ dunn.test(all_data$lesion,all_data$group, method="bonferroni")  # non-normalized sobol indices


#~ args=c('/media/simulations/WML/analysis/CSV/M1ROI/FAZEKAS1_M1_ROI_left_stats_summary.csv','/media/simulations/WML/analysis/CSV/M1ROI/FAZEKAS2_M1_ROI_left_stats_summary.csv','/media/simulations/WML/analysis/CSV/M1ROI/FAZEKAS3_M1_ROI_left_stats_summary.csv')
#~ reader<-DataProvider(args,TRUE)
#~ all_data<-reader$stacked_input
#~ all_data$subjectID<-NULL
#~ pairwise.t.test(all_data$scalp,all_data$group, p.adj="bonf")     # normalized sobol indices
#~ pairwise.t.test(all_data$gm,all_data$group, p.adj="bonf")        # normalized sobol indices
#~ dunn.test(all_data$wm,all_data$group, method="bonferroni")       # normalized sobol indices
#~ dunn.test(all_data$lesion,all_data$group, method="bonferroni")   # normalized sobol indices




#~ args=c('/media/simulations/WML/analysis/CSV/M1ROI/FAZEKAS1_M1_ROI_right_stats_summary.csv','/media/simulations/WML/analysis/CSV/M1ROI/FAZEKAS2_M1_ROI_right_stats_summary.csv','/media/simulations/WML/analysis/CSV/M1ROI/FAZEKAS3_M1_ROI_right_stats_summary.csv')
# args=c('/media/simulations/WML/analysis/CSV/M1ROI_atlasbased/FAZEKAS1_M1_atlasbased_right_stats_summary.csv','/media/simulations/WML/analysis/CSV/M1ROI_atlasbased/FAZEKAS2_M1_atlasbased_right_stats_summary.csv','/media/simulations/WML/analysis/CSV/M1ROI_atlasbased/FAZEKAS3_M1_atlasbased_right_stats_summary.csv')
args=c('/media/simulations/WML/analysis/CSV/non-normalized_sobol_indices/FAZEKAS0_M1_atlasbased_right_stats_summary.csv','/media/simulations/WML/analysis/CSV/non-normalized_sobol_indices/FAZEKAS1_M1_atlasbased_right_stats_summary.csv','/media/simulations/WML/analysis/CSV/non-normalized_sobol_indices/FAZEKAS2_M1_atlasbased_right_stats_summary.csv','/media/simulations/WML/analysis/CSV/non-normalized_sobol_indices/FAZEKAS3_M1_atlasbased_right_stats_summary.csv')
reader<-DataProvider(args,TRUE)
all_data<-reader$stacked_input
all_data$subjectID<-NULL
#~ pairwise.t.test(all_data$scalp,all_data$group, p.adj="bonf")        # normalized sobol indices
#print("M1 right: wm")
#dunn.test(all_data$wm,all_data$group, method="bonferroni")          # normalized sobol indices
#print("M1 right: lesion")
#dunn.test(all_data$lesion,all_data$group, method="bonferroni")      # normalized sobol indices
print("M1 right: wm")
dunn.test(all_data$wm,all_data$group, method="bonferroni")          # non-normalized sobol indices
print("M1 right: lesion")
dunn.test(all_data$lesion,all_data$group, method="bonferroni")      # non-normalized sobol indices


#~ args=c('/media/simulations/WML/analysis/CSV/M1ROI/FAZEKAS1_M1_ROI_left_stats_summary.csv','/media/simulations/WML/analysis/CSV/M1ROI/FAZEKAS2_M1_ROI_left_stats_summary.csv','/media/simulations/WML/analysis/CSV/M1ROI/FAZEKAS3_M1_ROI_left_stats_summary.csv')
#args=c('/media/simulations/WML/analysis/CSV/M1ROI_atlasbased/FAZEKAS1_M1_atlasbased_left_stats_summary.csv','/media/simulations/WML/analysis/CSV/M1ROI_atlasbased/FAZEKAS2_M1_atlasbased_left_stats_summary.csv','/media/simulations/WML/analysis/CSV/M1ROI_atlasbased/FAZEKAS3_M1_atlasbased_left_stats_summary.csv')
args=c('/media/simulations/WML/analysis/CSV/non-normalized_sobol_indices/FAZEKAS0_M1_atlasbased_left_stats_summary.csv','/media/simulations/WML/analysis/CSV/non-normalized_sobol_indices/FAZEKAS1_M1_atlasbased_left_stats_summary.csv','/media/simulations/WML/analysis/CSV/non-normalized_sobol_indices/FAZEKAS2_M1_atlasbased_left_stats_summary.csv','/media/simulations/WML/analysis/CSV/non-normalized_sobol_indices/FAZEKAS3_M1_atlasbased_left_stats_summary.csv')
reader<-DataProvider(args,TRUE)
all_data<-reader$stacked_input
all_data$subjectID<-NULL
#~ pairwise.t.test(all_data$scalp,all_data$group, p.adj="bonf")     # normalized sobol indices
#print("M1 left: wm")
#dunn.test(all_data$wm,all_data$group, method="bonferroni")       # normalized sobol indices
#print("M1 left: lesion")
#dunn.test(all_data$lesion,all_data$group, method="bonferroni")   # normalized sobol indices
#~ print("M1 left: scalp")
#~ pairwise.t.test(all_data$scalp,all_data$group, p.adj="bonf")        # non-normalized sobol indices
#~ print("M1 left: skull")
#~ pairwise.t.test(all_data$skull,all_data$group, p.adj="bonf")        # non-normalized sobol indices
print("M1 left: wm")
dunn.test(all_data$wm,all_data$group, method="bonferroni")          # non-normalized sobol indices
print("M1 left: lesion")
dunn.test(all_data$lesion,all_data$group, method="bonferroni")      # non-normalized sobol indices

#~ args=c('/media/simulations/WML/analysis/CSV/ElectrodeROI/FAZEKAS1_ROI_1_stats_summary.csv','/media/simulations/WML/analysis/CSV/ElectrodeROI/FAZEKAS2_ROI_1_stats_summary.csv','/media/simulations/WML/analysis/CSV/ElectrodeROI/FAZEKAS3_ROI_1_stats_summary.csv')
args=c('/media/simulations/WML/analysis/CSV/non-normalized_sobol_indices/FAZEKAS0_ROI1_stats_summary.csv','/media/simulations/WML/analysis/CSV/non-normalized_sobol_indices/FAZEKAS1_ROI1_stats_summary.csv','/media/simulations/WML/analysis/CSV/non-normalized_sobol_indices/FAZEKAS2_ROI1_stats_summary.csv','/media/simulations/WML/analysis/CSV/non-normalized_sobol_indices/FAZEKAS3_ROI1_stats_summary.csv')
reader<-DataProvider(args,TRUE)
all_data<-reader$stacked_input
all_data$subjectID<-NULL
#print("ROI1: scalp")
#pairwise.t.test(all_data$scalp,all_data$group, p.adj="bonf")    # normalized sobol indices
#print("ROI1: wm")
#dunn.test(all_data$wm,all_data$group, method="bonferroni")      # normalized sobol indices
#print("ROI1: lesion")
#dunn.test(all_data$lesion,all_data$group, method="bonferroni")  # normalized sobol indices

print("ROI1: scalp")
dunn.test(all_data$scalp,all_data$group, method="bonferroni")   # non-normalized sobol indices
print("ROI1: wm")
dunn.test(all_data$wm,all_data$group, method="bonferroni")      # non-normalized sobol indices
print("ROI1: lesion")
dunn.test(all_data$lesion,all_data$group, method="bonferroni")  # non-normalized sobol indices


#args=c('/media/simulations/WML/analysis/CSV/ElectrodeROI/FAZEKAS1_ROI_2_stats_summary.csv','/media/simulations/WML/analysis/CSV/ElectrodeROI/FAZEKAS2_ROI_2_stats_summary.csv','/media/simulations/WML/analysis/CSV/ElectrodeROI/FAZEKAS3_ROI_2_stats_summary.csv')
args=c('/media/simulations/WML/analysis/CSV/non-normalized_sobol_indices/FAZEKAS0_ROI2_stats_summary.csv','/media/simulations/WML/analysis/CSV/non-normalized_sobol_indices/FAZEKAS1_ROI2_stats_summary.csv','/media/simulations/WML/analysis/CSV/non-normalized_sobol_indices/FAZEKAS2_ROI2_stats_summary.csv','/media/simulations/WML/analysis/CSV/non-normalized_sobol_indices/FAZEKAS3_ROI2_stats_summary.csv')
reader<-DataProvider(args,TRUE)
all_data<-reader$stacked_input
all_data$subjectID<-NULL
#print("ROI2: scalp")
#pairwise.t.test(all_data$scalp,all_data$group, p.adj="bonf")    # normalized sobol indices
#print("ROI2: wm")
#dunn.test(all_data$wm,all_data$group, method="bonferroni")      # normalized sobol indices
#print("ROI2: lesion")
#dunn.test(all_data$lesion,all_data$group, method="bonferroni")  # normalized sobol indices

print("ROI2: wm")
dunn.test(all_data$wm,all_data$group, method="bonferroni")      # non-normalized sobol indices
print("ROI2: lesion")
dunn.test(all_data$lesion,all_data$group, method="bonferroni")  # non-normalized sobol indices
#~ print("ROI2: wm")
#~ dunn.test(all_data$wm,all_data$group, method="bonferroni")      # non-normalized sobol indices
