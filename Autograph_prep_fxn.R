###this is a script to make the analysis of FlowJo data easier, basically it will take the .xls file
### that is output by FlowJo and turn it  into something the FC_statistics autograph fxn can use

#the first argument in this function is sample_table, and it's just what it sounds like: a description
# of the samples, such as age and whether they are WT or Tg. The Sample name MUST be the same as the
# one in FlowJo. Must be in the format:

# Sample    Age    Genotype
#  B1       2       Tg
#  B2       12      WT

AutoGraph_prep <- function(file_path, sample_table_path, CD45_Phrase, CD45Hi_Phrase) {
  library(tidyr)
  library(dplyr)

  #load in the sample table
  sample_table <- read.delim(file = sample_table_path, header = TRUE, sep = "\t")

  #then let's load in the data
  data <- read.delim(file = file_path, header = TRUE, sep = "\t")
  #next, let's split the name column at the first period to produce a sample names column:
  data <- extract(data, Name, into = c("Name", "X"), "^([^.]+)\\.(.*)")

  #now, lets add stuff to the sample_Table that we want to include in adding to the file
  cd45_df <- subset(data, endsWith(data$X, CD45_Phrase))
  sample_table$CD45_number <- cd45_df$Total.Cell.Number
  cd45hi_df <- subset(data, endsWith(data$X, CD45Hi_Phrase))
  sample_table$CD45Hi_number <- cd45hi_df$Total.Cell.Number

  #now let's add all this information to the data object
  data <- left_join(data, sample_table, by = c("Name" = "Sample"))

  #Now we need to normalize the total cell number by the CD45 column
  data$CD45_norm <- data$Total.Cell.Number / data$CD45_number * 100
  #next normalize the total cell number by the CD45 hi column
  data$CD45Hi_norm <- data$Total.Cell.Number / data$CD45Hi_number * 100

  #get rid of the extraneous columns
  data$CD45_number <- NULL
  data$CD45Hi_number <- NULL
  #Rename the column to something that makes sense
  names(data)[4] <- "Percentage.of.Parent"

  #Get rid of any numbers over 100 that were created by normalization
  data[8:9] <- lapply(data[8:9], function(x) ifelse(x>100.1, NA, x))

  #reorganize the dataframe into the desired format
  data <- data[, c("Genotype", "Age", "Name", "X", "Percentage.of.Parent", "Total.Cell.Number",
                   "CD45_norm", "CD45Hi_norm")]
  return(data)
}
