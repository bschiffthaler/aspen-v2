#' ---
#' title: "Updating the ENA sample information to add new mandatory fields"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(xml2)
})

#' # Input
metadata <- read_tsv(here("doc/ENA/updated_biosample_metadata_Robinson.txt"),
                     name_repair = "universal",show_col_types = FALSE) %>% 
  select(Sample.ID,latitudeN,longitudeE,collection.date) %>% 
  mutate(collection.date=as.character(collection.date))

xml <- read_xml(here("doc/ENA/UPSC-0182.Sample.xml"))

#' # Process
children <- xml_children(xml)

dev.null <- lapply(1:length(children),function(i,ch){

  p <- which(metadata$Sample.ID == children[[i]] %>% xml_child() %>% xml_text())
  
  children[[i]] %>% xml_add_child("SAMPLE_ATTRIBUTES") %>% 
    xml_add_child("SAMPLE_ATTRIBUTE") %>% 
    xml_add_child("TAG","geographic location (country and/or sea)") %>% 
    xml_add_sibling("VALUE","Sweden") %>% 
    xml_parent() %>% xml_parent() %>% 
    xml_add_child("SAMPLE_ATTRIBUTE") %>% 
    xml_add_child("TAG","collection date") %>% 
    xml_add_sibling("VALUE",metadata$collection.date[[p]]) %>% 
    xml_parent() %>% xml_parent() %>% 
    xml_add_child("SAMPLE_ATTRIBUTE") %>% 
    xml_add_child("TAG","geographic location (latitude)") %>% 
    xml_add_sibling("VALUE",metadata$latitudeN[[p]]) %>% 
    xml_add_sibling("UNITS","DD") %>% 
    xml_parent() %>% xml_parent() %>% 
    xml_add_child("SAMPLE_ATTRIBUTE") %>% 
    xml_add_child("TAG","geographic location (longitude)") %>% 
    xml_add_sibling("VALUE",metadata$longitudeE[[p]]) %>% 
    xml_add_sibling("UNITS","DD")
},children)

#' # Validate
#' This is a purely visual check
children[[1]] %>% xml_child("SAMPLE_ATTRIBUTES")
xml_structure(children[[1]])
xml_structure(children[[242]])

#' # Export
write_xml(children %>% xml_root(),"~/OneDrive - Sveriges lantbruksuniversitet/UPSC/Facility/Doc/Tiggy/UPSC-0182.Updated.Sample.xml")

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

