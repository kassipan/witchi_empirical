#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#Left join to join the table by selecting all the records from the first dataframe and only matching records in the second dataframe.

in_table1 <- args[1]
in_table2 <- args[2]
joinByColumn <- args[3]
outputFile <- args[4]

# load packages
library(data.table)
library(dplyr)

x_table <- fread(in_table1, header = T, sep = "\t", fill = T)

y_table <- fread(in_table2, header = T, sep = "\t", fill = T)

left_joined_table <- left_join(x_table, y_table, by = joinByColumn)


write.table(left_joined_table, file = outputFile, row.names=FALSE, quote = FALSE, sep="\t")

#Example run
#Rscript --vanilla /local/two/panag007/scripts/left_join_files.R x_table y_table Bin_Id left_joined_table
