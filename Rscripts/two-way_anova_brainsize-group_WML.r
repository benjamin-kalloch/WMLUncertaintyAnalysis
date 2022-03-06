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
#library(ARTool)
library(hash)

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
                input[[i]] <- tryCatch(read.csv( file=input_CSVs[[i]] ), error=function(e){print("File does not exist! Skipping input file"); print(input_CSVs[[i]]);print(e)})

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

input_fnames=hash()
input_fnames[["oz-fpz_wholebrain"]]<-c('H:/PyGPC_sims/WML/analysis_oz_fpz/CSV_summary_non-normalized_sobol_indices_with_brain_sizes/FAZEKAS0_wholebrain_stats_summary_oz_fpz.csv','H:/PyGPC_sims/WML/analysis_oz_fpz/CSV_summary_non-normalized_sobol_indices_with_brain_sizes/FAZEKAS1_wholebrain_stats_summary_oz_fpz.csv','H:/PyGPC_sims/WML/analysis_oz_fpz/CSV_summary_non-normalized_sobol_indices_with_brain_sizes/FAZEKAS2_wholebrain_stats_summary_oz_fpz.csv','H:/PyGPC_sims/WML/analysis_oz_fpz/CSV_summary_non-normalized_sobol_indices_with_brain_sizes/FAZEKAS3_wholebrain_stats_summary_oz_fpz.csv')
input_fnames[["oz-fpz_electrode_ROI_right"]]<-c('H:/PyGPC_sims/WML/analysis_oz_fpz/CSV_summary_non-normalized_sobol_indices_with_brain_sizes/FAZEKAS0_ROI_2_stats_summary_oz_fpz.csv','H:/PyGPC_sims/WML/analysis_oz_fpz/CSV_summary_non-normalized_sobol_indices_with_brain_sizes/FAZEKAS1_ROI_2_stats_summary_oz_fpz.csv','H:/PyGPC_sims/WML/analysis_oz_fpz/CSV_summary_non-normalized_sobol_indices_with_brain_sizes/FAZEKAS2_ROI_2_stats_summary_oz_fpz.csv','H:/PyGPC_sims/WML/analysis_oz_fpz/CSV_summary_non-normalized_sobol_indices_with_brain_sizes/FAZEKAS3_ROI_2_stats_summary_oz_fpz.csv')
input_fnames[["oz-fpz_electrode_ROI_left"]]<-c('H:/PyGPC_sims/WML/analysis_oz_fpz/CSV_summary_non-normalized_sobol_indices_with_brain_sizes/FAZEKAS0_ROI_1_stats_summary_oz_fpz.csv','H:/PyGPC_sims/WML/analysis_oz_fpz/CSV_summary_non-normalized_sobol_indices_with_brain_sizes/FAZEKAS1_ROI_1_stats_summary_oz_fpz.csv','H:/PyGPC_sims/WML/analysis_oz_fpz/CSV_summary_non-normalized_sobol_indices_with_brain_sizes/FAZEKAS2_ROI_1_stats_summary_oz_fpz.csv','H:/PyGPC_sims/WML/analysis_oz_fpz/CSV_summary_non-normalized_sobol_indices_with_brain_sizes/FAZEKAS3_ROI_1_stats_summary_oz_fpz.csv')
input_fnames[["oz-fpz_M1_ROI_right"]]<-c('H:/PyGPC_sims/WML/analysis_oz_fpz/CSV_summary_non-normalized_sobol_indices_with_brain_sizes/FAZEKAS0_M1_ROI_right_stats_summary_oz_fpz.csv','H:/PyGPC_sims/WML/analysis_oz_fpz/CSV_summary_non-normalized_sobol_indices_with_brain_sizes/FAZEKAS1_M1_ROI_right_stats_summary_oz_fpz.csv','H:/PyGPC_sims/WML/analysis_oz_fpz/CSV_summary_non-normalized_sobol_indices_with_brain_sizes/FAZEKAS2_M1_ROI_right_stats_summary_oz_fpz.csv','H:/PyGPC_sims/WML/analysis_oz_fpz/CSV_summary_non-normalized_sobol_indices_with_brain_sizes/FAZEKAS3_M1_right_stats_summary_oz_fpz.csv')
input_fnames[["oz-fpz_M1_ROI_left"]]<-c('H:/PyGPC_sims/WML/analysis_oz_fpz/CSV_summary_non-normalized_sobol_indices_with_brain_sizes/FAZEKAS0_M1_ROI_left_stats_summary_oz_fpz.csv','H:/PyGPC_sims/WML/analysis_oz_fpz/CSV_summary_non-normalized_sobol_indices_with_brain_sizes/FAZEKAS1_M1_ROI_left_stats_summary_oz_fpz.csv','H:/PyGPC_sims/WML/analysis_oz_fpz/CSV_summary_non-normalized_sobol_indices_with_brain_sizes/FAZEKAS2_M1_ROI_left_stats_summary_oz_fpz.csv','H:/PyGPC_sims/WML/analysis_oz_fpz/CSV_summary_non-normalized_sobol_indices_with_brain_sizes/FAZEKAS3_M1_left_stats_summary_oz_fpz.csv')
input_fnames[["oz-fpz_hippocampus_L"]]<-c('H:/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/oz-fpz/with_brainsize/FAZEKAS0_hippocampus_L_stats_summary_oz-fpz.csv','H:/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/oz-fpz/with_brainsize/FAZEKAS1_hippocampus_L_stats_summary_oz-fpz.csv','H:/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/oz-fpz/with_brainsize/FAZEKAS2_thalamus_L_stats_summary_oz-fpz.csv','H:/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/oz-fpz/with_brainsize/FAZEKAS3_hippocampus_L_stats_summary_oz-fpz.csv')
input_fnames[["oz-fpz_hippocampus_R"]]<-c('H:/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/oz-fpz/with_brainsize/FAZEKAS0_hippocampus_R_stats_summary_oz-fpz.csv','H:/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/oz-fpz/with_brainsize/FAZEKAS1_hippocampus_R_stats_summary_oz-fpz.csv','H:/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/oz-fpz/with_brainsize/FAZEKAS2_thalamus_R_stats_summary_oz-fpz.csv','H:/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/oz-fpz/with_brainsize/FAZEKAS3_hippocampus_R_stats_summary_oz-fpz.csv')
input_fnames[["oz-fpz_thalamus_L"]]<-c('H:/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/oz-fpz/with_brainsize/FAZEKAS0_thalamus_L_stats_summary_oz-fpz.csv','H:/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/oz-fpz/with_brainsize/FAZEKAS1_thalamus_L_stats_summary_oz-fpz.csv','H:/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/oz-fpz/with_brainsize/FAZEKAS2_thalamus_L_stats_summary_oz-fpz.csv','H:/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/oz-fpz/with_brainsize/FAZEKAS3_thalamus_L_stats_summary_oz-fpz.csv')
input_fnames[["oz-fpz_thalamus_R"]]<-c('H:/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/oz-fpz/with_brainsize/FAZEKAS0_thalamus_R_stats_summary_oz-fpz.csv','H:/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/oz-fpz/with_brainsize/FAZEKAS1_thalamus_R_stats_summary_oz-fpz.csv','H:/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/oz-fpz/with_brainsize/FAZEKAS2_thalamus_R_stats_summary_oz-fpz.csv','H:/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/oz-fpz/with_brainsize/FAZEKAS3_thalamus_R_stats_summary_oz-fpz.csv')


#input_fnames[["c3-c4_wholebrain"]]<-c('/media/benny/work_backup/Promotion/PyGPC_sims/WML/analysis/CSV/non-normalized_sobol_indices/with_brainsize/FAZEKAS0_wholebrain_stats_summary.csv','/media/benny/work_backup/Promotion/PyGPC_sims/WML/analysis/CSV/non-normalized_sobol_indices/with_brainsize/FAZEKAS1_wholebrain_stats_summary.csv','/media/benny/work_backup/Promotion/PyGPC_sims/WML/analysis/CSV/non-normalized_sobol_indices/with_brainsize/FAZEKAS2_wholebrain_stats_summary.csv','/media/benny/work_backup/Promotion/PyGPC_sims/WML/analysis/CSV/non-normalized_sobol_indices/with_brainsize/FAZEKAS3_wholebrain_stats_summary.csv')
#input_fnames[["c3-c4_electrode_ROI_right"]]<-c('/media/benny/work_backup/Promotion/PyGPC_sims/WML/analysis/CSV/non-normalized_sobol_indices/with_brainsize/FAZEKAS0_ROI_2_stats_summary.csv','/media/benny/work_backup/Promotion/PyGPC_sims/WML/analysis/CSV/non-normalized_sobol_indices/with_brainsize/FAZEKAS1_ROI2_stats_summary.csv','/media/benny/work_backup/Promotion/PyGPC_sims/WML/analysis/CSV/non-normalized_sobol_indices/with_brainsize/FAZEKAS2_ROI2_stats_summary.csv','/media/benny/work_backup/Promotion/PyGPC_sims/WML/analysis/CSV/non-normalized_sobol_indices/with_brainsize/FAZEKAS3_ROI2_stats_summary.csv')
#input_fnames[["c3-c4_electrode_ROI_left"]]<-c('/media/benny/work_backup/Promotion/PyGPC_sims/WML/analysis/CSV/non-normalized_sobol_indices/with_brainsize/FAZEKAS0_ROI_1_stats_summary.csv','/media/benny/work_backup/Promotion/PyGPC_sims/WML/analysis/CSV/non-normalized_sobol_indices/with_brainsize/FAZEKAS1_ROI1_stats_summary.csv','/media/benny/work_backup/Promotion/PyGPC_sims/WML/analysis/CSV/non-normalized_sobol_indices/with_brainsize/FAZEKAS2_ROI1_stats_summary.csv','/media/benny/work_backup/Promotion/PyGPC_sims/WML/analysis/CSV/non-normalized_sobol_indices/with_brainsize/FAZEKAS3_ROI1_stats_summary.csv')
#input_fnames[["c3-c4_M1_ROI_right"]]<-c('/media/benny/work_backup/Promotion/PyGPC_sims/WML/analysis/CSV/non-normalized_sobol_indices/with_brainsize/FAZEKAS0_M1_right_stats_summary.csv','/media/benny/work_backup/Promotion/PyGPC_sims/WML/analysis/CSV/non-normalized_sobol_indices/with_brainsize/FAZEKAS1_M1_right_stats_summary.csv','/media/benny/work_backup/Promotion/PyGPC_sims/WML/analysis/CSV/non-normalized_sobol_indices/with_brainsize/FAZEKAS2_M1_right_stats_summary.csv','/media/benny/work_backup/Promotion/PyGPC_sims/WML/analysis/CSV/non-normalized_sobol_indices/with_brainsize/FAZEKAS3_M1_right_stats_summary.csv')
#input_fnames[["c3-c4_M1_ROI_left"]]<-c('/media/benny/work_backup/Promotion/PyGPC_sims/WML/analysis/CSV/non-normalized_sobol_indices/with_brainsize/FAZEKAS0_M1_left_stats_summary.csv','/media/benny/work_backup/Promotion/PyGPC_sims/WML/analysis/CSV/non-normalized_sobol_indices/with_brainsize/FAZEKAS1_M1_left_stats_summary.csv','/media/benny/work_backup/Promotion/PyGPC_sims/WML/analysis/CSV/non-normalized_sobol_indices/with_brainsize/FAZEKAS2_M1_left_stats_summary.csv','/media/benny/work_backup/Promotion/PyGPC_sims/WML/analysis/CSV/non-normalized_sobol_indices/with_brainsize/FAZEKAS3_M1_left_stats_summary.csv')
#input_fnames[["c3-c4_hippocampus_L"]]<-c('/home/benny/NAS/htwk/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/c3-c4/with_brainsize/FAZEKAS0_hippocampus_L_stats_summary_c3-c4.csv','/home/benny/NAS/htwk/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/c3-c4/with_brainsize/FAZEKAS1_hippocampus_L_stats_summary_c3-c4.csv','/home/benny/NAS/htwk/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/c3-c4/with_brainsize/FAZEKAS2_thalamus_L_stats_summary_c3-c4.csv','/home/benny/NAS/htwk/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/c3-c4/with_brainsize/FAZEKAS3_hippocampus_L_stats_summary_c3-c4.csv')
#input_fnames[["c3-c4_hippocampus_R"]]<-c('/home/benny/NAS/htwk/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/c3-c4/with_brainsize/FAZEKAS0_hippocampus_R_stats_summary_c3-c4.csv','/home/benny/NAS/htwk/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/c3-c4/with_brainsize/FAZEKAS1_hippocampus_R_stats_summary_c3-c4.csv','/home/benny/NAS/htwk/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/c3-c4/with_brainsize/FAZEKAS2_thalamus_R_stats_summary_c3-c4.csv','/home/benny/NAS/htwk/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/c3-c4/with_brainsize/FAZEKAS3_hippocampus_R_stats_summary_c3-c4.csv')
#input_fnames[["c3-c4_thalamus_L"]]<-c('/home/benny/NAS/htwk/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/c3-c4/with_brainsize/FAZEKAS0_thalamus_L_stats_summary_c3-c4.csv','/home/benny/NAS/htwk/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/c3-c4/with_brainsize/FAZEKAS1_thalamus_L_stats_summary_c3-c4.csv','/home/benny/NAS/htwk/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/c3-c4/with_brainsize/FAZEKAS2_thalamus_L_stats_summary_c3-c4.csv','/home/benny/NAS/htwk/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/c3-c4/with_brainsize/FAZEKAS3_thalamus_L_stats_summary_c3-c4.csv')
#input_fnames[["c3-c4_thalamus_R"]]<-c('/home/benny/NAS/htwk/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/c3-c4/with_brainsize/FAZEKAS0_thalamus_R_stats_summary_c3-c4.csv','/home/benny/NAS/htwk/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/c3-c4/with_brainsize/FAZEKAS1_thalamus_R_stats_summary_c3-c4.csv','/home/benny/NAS/htwk/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/c3-c4/with_brainsize/FAZEKAS2_thalamus_R_stats_summary_c3-c4.csv','/home/benny/NAS/htwk/htwk/data/medical/LIFE/subjects/subcortical_roi/analysis_csv/c3-c4/with_brainsize/FAZEKAS3_thalamus_R_stats_summary_c3-c4.csv')


for(v in ls(input_fnames)){
    fnames<-input_fnames[[v]]
    print(v)
    reader <- DataProvider(fnames, TRUE)
    reader$stacked_input$subjectID <- NULL
    #reader$stacked_input <- reader$stacked_input[ ! reader$stacked_input$group == 0, ]
    reader$trimStackedColumns(c('lesion','group','brainsize', 'mean'))
    reader$trimmed_input$lesion<-reader$trimmed_input$lesion * 1e8

    reader$trimmed_input$brainsize_grouped <- cut(reader$trimmed_input$brainsize, breaks=5, dig.lab=10)
    print(table(reader$trimmed_input$brainsize_grouped))
        
    print("Kruskall-Wallis, one-way")
    reader$trimmed_input$brainsize_grouped<-as.numeric(reader$trimmed_input$brainsize_grouped)
    shap_test_res<-shapiro.test( reader$trimmed_input$brainsize ) # normality assumption violated if p < 0.05
    print(shap_test_res)
    print(cor.test(reader$trimmed_input$lesion, reader$trimmed_input$brainsize, method="kendall"))
    print(kruskal.test( mean ~ brainsize_grouped, reader$trimmed_input))
    
    
    
#    print("ANOVA")
#    two.way <- aov(lesion ~ group + brainsize, data=reader$trimmed_input)
#    interaction <- aov(lesion ~ group * brainsize, data=reader$trimmed_input)
#    two.way.interaction <- aov(lesion ~ group + brainsize + group * brainsize, data=reader$trimmed_input)
#    print(summary(two.way))
#    print(summary(interaction))
#    print(summary(two.way.interaction))
    
    
#    print("Aligned Rank Transform for Nonparametric Factorial ANOVAs")
#    two_way_with_interaction<-art(lesion ~ group + brainsize + group:brainsize, data=reader$trimmed_input)
#    two_way_with_interaction_aov<-anova(two_way_with_interaction)
#    print(summary(two_way_with_interaction))
#    print(two_way_with_interaction_aov)
    cat('\n\n\n\n')
}
