library(stringr)

#
# Usage: provide the paths of the CSv files containig the ROI data
#   - each file must contain all the data for this specfic subject.
#      Expected structure:
#       column 1: mean E-field value (magnitude)
#       column 2: variance of E-Field
#       column 3 - 23: sobol indices
#       column 24: point ID (ignored/skipped)
#   - each file is considered a unique subject.
#     The id of the subject must be part of the file name.
#     You can change the variable "position_of_subjectID_in_path"
#     to indicate where the ID can be found.
#
# Possibly need modification:
# - The regex identifying the sobol indices in 'strip_prefix_from_column_names'
# - The position of the subjectID in the path string
# - The name of the column containing the variance and the mean (if sobol indices are normalized and the normalization should be reversed)
# - The name of the output file
# - Note: We assume normalized Sobol indices in the raw data here!
#

# function definitions
strip_prefix_from_column_names<-function( df ){
    column_names <- colnames( df )

#    column_names <- str_match( column_names, 'ROI_left_(E_)?(sobol_[0-9]+_)?(.*)' ); matching_group = 4
#     column_names <- str_match( column_names, 'ROI_right_(E_)?(sobol_[0-9]+_)?(.*)' ); matching_group = 4
#~     column_names <- str_match( column_names, 'ROI_inner_(E_)?(sobol_[0-9]+_)?(.*)' ); matching_group = 4
     column_names <- str_match( column_names, 'E_(sobol_[0-9]+_)?(.*)' ); matching_group = 3

    colnames( df ) <- column_names[,matching_group]

    return( df )
}

remove_na_columns<-function( df ){
    df <- df[!is.na(names(df))]

    return(df)
}

create_row_for_result_dataframe<-function( data_frame, subjectID ) {
    data_frame <- as.data.frame( t(colMeans( data_frame )) )
    data_frame$subjectID <- subjectID

    return( data_frame )
}

multiply_sobol_indices_by_variance<-function( data_frame ) {
    column_names = names(data_frame)

    for( c in 1:ncol( data_frame ) ) {
         if( column_names[c] != "mean"  & column_names[c] != "var" ) {
             data_frame[c] = data_frame[c] * data_frame$var
         }
    }

    return( data_frame )
}

#### Main ####

# settings
POSITION_OF_SUBJECTID_IN_PATH=8
USE_NORMALIZED_SOBOL_INDICES=FALSE

# parse commandline arguments
args<-commandArgs(trailingOnly=TRUE)
if( length(args) < 1 ) {
    stop("Please provide the paths to the CSV files of all subjects, each containing the values inside the ROI per subject.", call.=FALSE)
}

# read the files in
roi_data<-lapply( args, read.csv, header=TRUE )


#header = c("subjectID","mean","var","gm","scalp","wm","skull","csf","scalp_gm","csf_gm","scalp_skull","scalp_csf","scalp_wm","lesion","csf_wm","wm_lesion","gm_wm","skull_wm","skull_gm","gm_lesion","csf_lesion","skull_csf","skull_lesion","scalp_lesion");
header = c("subjectID","mean","var","scalp","skull","csf","gm","wm","scalp_gm","csf_gm","scalp_skull","scalp_csf","scalp_wm","csf_wm","gm_wm","skull_wm","skull_gm","skull_csf");

subjects = vector( mode = "list", length=length(args) )
df = data.frame(matrix(ncol = length(header), nrow = 0))

# read the subject IDs from the input paths
for( i in 1:length(args) ) {
    subjects[[i]] = strsplit( args[[i]], "/" )[[1]][[POSITION_OF_SUBJECTID_IN_PATH]]
}
colnames( df ) <- header

for( i in 1:length(args) ) {
    roi_data[[i]] <- strip_prefix_from_column_names( roi_data[[i]] )
    roi_data[[i]] <- remove_na_columns( roi_data[[i]] )
    if( ! USE_NORMALIZED_SOBOL_INDICES ) {
        roi_data[[i]] <- multiply_sobol_indices_by_variance( roi_data[[i]] )
    }
    df <- rbind( df, create_row_for_result_dataframe( roi_data[[i]], subjects[[i]] ) )
}

write.csv( df, file = "/home/benny/ram_folder/FAZEKAS0_wholebrain_stats_summary_non-normalized_sobol.csv", row.names = TRUE )

print( df )

#print("Statistics over all ROIs (column_name:mean:variance:std_dev)")
#full_roi_data=Reduce(function(x,y){merge(x,y,all=TRUE)},roi_data)
#print_col_stats( full_roi_data )
