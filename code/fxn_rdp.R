# functions to format RDP results

# read in an RDP result table and tidy up the formatting -- output from the web interface
read.format.rdp_fromWeb <- function(filepath){
  
  rdp.result <- read.delim(filepath, header = 0, skip = 6, sep = ";", stringsAsFactors = F)
  new.colnames <- c("read", "dir",
                    "domain","d.perc",
                    "phylum","p.perc",
                    "class","c.perc",
                    "order","o.perc",
                    "family","f.perc",
                    "genus","g.perc",
                    "species", "s.perc")
  colnames(rdp.result) <- new.colnames
  
  return(rdp.result)
}

# read in an RDP result table and tidy up the formatting -- output from the command line interface
read.format.rdp_fromCL <- function(filepath){
  
  rdp.result <- read.delim(filepath, header = 0)
  removecols <- c(2, 4, 7, 10, 13, 16, 19, 22)
  new.colnames <- c("read",
                    "domain","d.perc",
                    "phylum","p.perc",
                    "class","c.perc",
                    "order","o.perc",
                    "family","f.perc",
                    "genus","g.perc",
                    "species", "s.perc")
  rdp.result <- rdp.result[,-removecols]
  colnames(rdp.result) <- new.colnames
  
  return(rdp.result)
}

unclassified_at <- function(rdp.df, perc.conf.cutoff){
  
  rdp.df %>%
    separate(d.perc, into = c("d.perc","drop")) %>%
    separate(p.perc, into = c("p.perc","drop")) %>%
    separate(c.perc, into = c("c.perc","drop")) %>%
    separate(o.perc, into = c("o.perc","drop")) %>%
    separate(f.perc, into = c("f.perc","drop")) %>%
    separate(g.perc, into = c("g.perc","drop")) %>%
    select(-c(drop, s.perc)) %>%
    mutate(d.perc = as.numeric(d.perc)) %>%
    mutate(p.perc = as.numeric(p.perc)) %>%
    mutate(c.perc = as.numeric(c.perc)) %>%
    mutate(o.perc = as.numeric(o.perc)) %>%
    mutate(f.perc = as.numeric(f.perc)) %>%
    mutate(g.perc = as.numeric(g.perc)) -> rdp.df.tmp
  
  # replace anything less than X% with unclassified
  rdp.df.tmp %>%
    mutate(domain = ifelse(d.perc < perc.conf.cutoff, "unclassified", as.character(domain))) %>%
    mutate(phylum = ifelse(p.perc < perc.conf.cutoff, "unclassified", as.character(phylum))) %>%
    mutate(class = ifelse(c.perc < perc.conf.cutoff, "unclassified", as.character(class))) %>%
    mutate(order = ifelse(o.perc < perc.conf.cutoff, "unclassified", as.character(order))) %>%
    mutate(family = ifelse(f.perc < perc.conf.cutoff, "unclassified", as.character(family))) %>%
    mutate(genus = ifelse(g.perc < perc.conf.cutoff, "unclassified", as.character(genus))) %>%
    mutate(species = as.character(species)) -> tmp
  
  # fill in "unclassified" down the tree
  for(i in 1:dim(tmp)[1]){
    
    if(tmp[i,"domain"] == "unclassified"){
      tmp[i,c("phylum","class","order","family","genus","species")] <- "unclassified"
    }
    
    if(tmp[i,"phylum"] == "unclassified"){
      tmp[i,c("class","order","family","genus","species")] <- "unclassified"
    }
    
    if(tmp[i,"class"] == "unclassified"){
      tmp[i,c("order","family","genus","species")] <- "unclassified"
    }
    
    if(tmp[i,"order"] == "unclassified"){
      tmp[i,c("family","genus","species")] <- "unclassified"
    }
    
    if(tmp[i,"family"] == "unclassified"){
      tmp[i,c("genus","species")] <- "unclassified"
    }
    
    if(tmp[i,"genus"] == "unclassified"){
      tmp[i,c("species")] <- "unclassified"
    }
    
  }
  
  return(tmp)
  
}
