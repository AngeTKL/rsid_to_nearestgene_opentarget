#' gene_id search
#'
#' Main function of the function - returns gene HGNC names using OpenTargets.
#'
#'nearest_rsid_to_gene_opentarget<-function(dat, Path, name_file)
#'dat : input dataset must contain a column 'rsid' with rsid . Ensure that colnames includes 'rsid'
#'for example the structure should be like this : 
#'dat
#'
          #'rsid               |  ......
          #'rs1257             |  ......
          #'rs2569             |  ......         
          #'rs12584            |  ......
          #'rs58454454845685   |  ......
#'Path: location where the output file will be saved
#'name_file: name of the output file 
#'
#'the function outputs a .csv file with two columns 
#'$Gene - for the nearest gene  & 
#'$rsid - for the corresponding  rsid .
#'eg of output: 
            #'rsid            \  Gene
            #'rs1257          \  AAD
            #'rs2569          \  XXS
            #'rs12584         \  RRRA
            #'rs58454454845685\  na
#'
#'!!!!!!!!!!!!!!!!!!!!!NOTE :  IF THe rsid in not in OPENTARGET database, the Gene will be "na"; !!!!!!!!!!!
#'by Ange Tchuisseu 
#'contact me if needed 
#'
  
nearest_gene_with_rsid_opentarget<-function(dat, Path, name_file)
  
  {

library(httr)
  message('print the number of rsid')
print (nrow(dat))
#check_required_columns(dat, "rsid")
#transform rsid column as a vector 
list_rsid <- as.vector(dat$rsid)
data= data.frame() #create an empty dataset

#######################################start the loop #################################
#loop to find the nearest gene name for each rsid based on opentarget
for (i in 1:nrow(dat)) {
     #nrow(list_rsid) ) {

# Set rs ID variable
query_rsID=list_rsid[i]
#query_rsID="rs114818383"

# Build query string TO FIND the nearest gene
query_string = "
query useSearchToConvertRSIDIntoIDFormat($query_rsID: String!) {
  search(queryString: $query_rsID) {
    variants {
      id
      rsId
      nearestGene { 
        id
        start
        symbol
        tss
        description
        chromosome
        exons
      }
      nearestGeneDistance
    }
  }
} 
"

# Set base URL of GraphQL API endpoint
base_url <- "https://api.genetics.opentargets.org/graphql"

# Set variables object of arguments to be passed to endpoint
variables <- list("query_rsID" = query_rsID)

# Construct POST request body object with query string and variables
post_body <- list(query = query_string, variables = variables)

# Perform POST request
r <- POST(url=base_url, body=post_body, encode='json')

# Results:
content(r)


df = content(r)

##########if the rsid is found in Opentarget then:

if (length(df$data$search$variants) !=0 ) { 
# Print first entry of nearest genes
head(df$data$search$variants[[1]]$nearestGene$symbol , 1)
# Flatten the nested result fields into a dataframe

list_result = as.data.frame(df$data$search$variants[[1]]$nearestGene$symbol)
colnames(list_result)='Gene' #rename the column
snp_id= as.data.frame(query_rsID)
colnames(snp_id)= 'rsid' #rename the column
x = cbind(list_result, snp_id ) #create a one line dataframe with one rsid and its corresponding nearest gene 
data=rbind(data, x[ 1,])  # merge with the other lines of rsid and genes 
}
##########if the rsid is not found in Opentarget then:
else { list_result = as.data.frame('na') #attribute the value 'na' as the corresponding gene
colnames(list_result)='Gene' #rename the column
snp_id= as.data.frame(query_rsID)
colnames(snp_id)= 'rsid' #rename the column
x = cbind(list_result, snp_id ) #create a one line dataframe with one rsid and its corresponding nearest gene 
data=rbind(data, x[ 1,]) # merge with the other lines of rsid and genes 
}

#the loop will continue up to nrow(dat)

}
#############################################end of the loop ##################
#to double check how many snps you have in the output file. 
message('The final dataset has the following number of lines', nrow(data))
print(nrow(data))
#to double check how many snps without corresponding gene are  in the output file ?
message('The final dataset has the following number of rsid not mapped', nrow(subset(data, data$Gene=='na')))
print(nrow(subset(data, data$Gene=='na')))

##############################output and save the result in a csv file ############################
write.csv(data, paste0(Path, '/', name_file, '.csv', sep=''), row.names = FALSE) 
}


#end 

