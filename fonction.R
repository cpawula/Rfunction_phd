#************************************************************************************************************************
## Script R
## Derniere modification : par Clovis Pawula
## Auteurs : 
#*************************************************************************************************************************#
#**************************************************** Fonction pour analyse données ssrseq ********************************#
#*************************************************************************************************************************#
##
##
##
##
##
#***********************************************************************************************************************


## Description du script :

# Ce script est utilisé pour rassembler les fonctions créées pour l'analyse des données ssrseq.




#***********************************************************************************************************************
# Dependencies:
# Package du tydiverse (dplyr, stringr, stringi etc..)
# Package polysat
# 
# Chargement des dependence


#***********************************************************************************************************************
#####Formatage des données#####

reformat_ssrseq_gt<-function(tab){
  # Cette fonction change le format du tableau de genotype donne en sortie du script de genotypage des individus
  tab_reform<-select(tab,1)%>%
    add_column(pop = "pls_enter", region = "pls_enter")
  for (i in seq(2,ncol(tab))) {
    newcolname<-str_replace_all(colnames(tab)[i],c("SSRseq_Rosa_" = "","_1._2._3._4" = ""))
    newcolnames<-c(paste0(newcolname,"_1"),paste0(newcolname,"_2"),paste0(newcolname,"_3"),paste0(newcolname,"_4"))
    splited<-str_split_fixed(tab[,i],"/",4)
    tab_reform[newcolnames]<-splited
  }
  return(tab_reform)
}


# Fonction qui, a partir, d'un format genind sort les deux fichiers necessaire à darwin le .don et le .var 
# En commentaire ci-dessous des lignes qui peuvent aider à passer d'un tableau en .txt a un objet genind utilisable dans la fonction
# tab<-read.table("GenotypicTable_AlleleDosage.txt", header = TRUE)%>%
#   as_tibble()%>%
#   mutate_if(is.integer, as.character)
# colnames(tab)<-str_remove_all(colnames(tab),"\\.")
# matrix<-df2genind(select(tab,-1),sep = "/",ploidy = 4, ind.names = tab$ind)
# matrix_to_Darwin(matrix,file = "opti_missing/darwin_format_for_opti")


matrix_to_Darwin<-function(genind, file, prsabs = FALSE){
  
  if (prsabs == TRUE) {
    genind@tab[genind@tab>1]<-1
  }
  
  var_RepA<-as.data.frame(tab(genind))%>%
    rownames_to_column()%>%
    dplyr::rename(Unit= rowname)
  var_RepA[is.na(var_RepA)] = 999
  
  
  Don_RepA<-as.data.frame(var_RepA$Unit)%>%
    rownames_to_column()
  var_RepA$Unit<-rownames(var_RepA)
  colnames(Don_RepA)<-c("Unit", "Geno")
  if (prsabs == TRUE) {
    line1var<-"@DARwin 5.0 - SINGLE"
  }else{line1var<-"@DARwin 5.0 - ALLELIC - 2"}
  
  line2var<-paste0(nrow(var_RepA),"\t",ncol(var_RepA))
  write(c(line1var,line2var), file = paste0(file,".VAR"),sep = "\t")
  write.table(var_RepA, file = paste0(file,".VAR"), row.names = FALSE,col.names = TRUE,quote = FALSE,sep = "\t",append = TRUE)
  
  line1don<-"@DARwin 5.0 - DON"
  write(c(line1don,line2var), file = paste0(file,".DON"),sep = "\t")
  write.table(Don_RepA, file = paste0(file,".DON"), row.names = FALSE,col.names = TRUE,quote = FALSE,sep = "\t",append = TRUE)
}



Quickboard<-function(data){
  write.table(data, "clipboard", sep="\t", col.names=TRUE, quote = F,row.names = FALSE)
}
# https://gist.github.com/zkamvar/aeaff83b9d126d55aade
gen2polysat <- function(gen, newploidy = gen@ploidy){
  if (!require(polysat)){
    stop("User needs polysat installed")
  }
  gen   <- recode_polyploids(gen, newploidy)
  gendf <- genind2df(gen, sep = "/", usepop = FALSE)
  gendf <- lapply(gendf, strsplit, "/")
  gendf <- lapply(gendf, lapply, as.numeric)
  ambig <- new("genambig", samples = indNames(gen), loci = locNames(gen))
  for (i in names(gendf)){
    res <- lapply(gendf[[i]], function(x) ifelse(is.na(x), Missing(ambig), x))
    Genotypes(ambig, loci = i) <- res
  }
  return(ambig)
}

Clones_removal<-function(genind, threshold = 0){
  # Cette fonction sert à retirer les clones pour ne garder qu'un individus représentatif 
  # l'individu retenu sera celui avc le moins de données manquantes, il est indiqué dans la dataframe ind_to_keep
  # une list des trois dataframe est retournée
  # info_clone = info des distances pairées
  # Control = info_clone mais avec uniquemetn les individus control
  # ind_to_keep = individus à conserver si on veut prendre celui avec le moins de missing data
  # Dependences : Polysat, Poppr, tidyverse
  # Modifée le 21 09 2022 par Clovis
  library(polysat)
  library(poppr)
  library(tidyverse)
  diagNA <- function(matrix) {
    
    diag(matrix)=NA
    return(matrix)
  }
  
  dist_ind<-dist(genind, method = "canberra", diag = TRUE)
  
  
  
  dist_ind<-as.matrix(dist_ind)
  clone_name<-assignClones(dist_ind,threshold = threshold)%>%as.data.frame()%>%
    rownames_to_column()%>%
    reshape::rename(c("rowname"="ind","."="clone" ))%>%
    rownames_to_column()%>%
    mutate(rowname = as.integer(rowname))
  
  # On calcule le taux de données manquantes de chaque individus pour ne conserver que l'individu ayant le moins de données manquantes par groupe de Clone
  clone_name$missing<-1-propTyped(genind,by="ind") 
  
  
  # On va produire un tableau récapitulant les distances entre paire
  a<-dist_ind%>%
    diagNA()
  info_clone<-as.data.frame(as.table(a))%>%
    dplyr::filter(!is.na(Freq))%>%
    mutate(comp = paste(Var2,Var1,sep = "_"))%>%
    mutate(Clone = Freq < threshold)
  
  control<-info_clone%>%
    dplyr::filter(str_detect(Var1,"_[AB]$") & str_detect(Var2,"_[AB]$"))%>%
    mutate(rep = str_match(Var1, str_sub(Var2,1,-2))[,1])
  
  # On produit un data.frame de noms d'individus à conserver 
  ind_to_keep<-clone_name%>%
    group_by(clone)%>%
    slice_min(missing,n = 1,with_ties = FALSE)
  
  output<-list(info_clone,control,ind_to_keep,clone_name)
  names(output) <- c("info_clone", "control","ind_to_keep", "clone_group")
  
  return(output)
}


recall_genind <- function(allele_tab,correspondance, prs_abs = FALSE) {
  a<-genind2df(allele_tab,sep = "_")
  i<-2
  # allele_tab : un objet genind
  # Un data.frame des correspondances entre code allélique
  # Cette fonction permet de changer le code allélique des individus 
  # Par exemple en remplacant les code des séquence par les code des tailles
  # Dependance : Poppr
  for (i in seq(2,length(colnames(a)))) {
    mark_unique<-filter(correspondance,str_detect(Locus,colnames(a)[i]))
    replacement<-mark_unique$AlleleLength
    names(replacement)<-mark_unique$AlleleSeqCode
    
    call_lgth<-str_replace_all(as.vector(a[,i]),replacement)
    a[,i]<-call_lgth
    i<-i+1
  }
  pop<-a$pop
  a<-a%>%
    dplyr::select(-pop)
  a<-df2genind(X = a,
               sep = "_",
               ploidy = allele_tab@ploidy,
               pop = allele_tab@pop,
               strata = allele_tab@strata)
  a@tab
  if (prs_abs == TRUE) {
    a@tab[a@tab>1]<-1
  }
  
  return(a)
}

VAR_to_genind<-function(VAR_tab,DON_tab){
  #Cette fonction sert à obtenir un objet genind à partir d'un .VAR et d'un .DON
  #Le .DON sert à obtenir le noms des individus
  #L'objet genind final n'a pas de pop associé ni de strata etc ...
  #VAR_tab : un data_frame issu d'un .VAR chargeable avec : read.table("RGal2012_tout20200401_CorrectionDM.VAR",header = T, sep = "\t", skip = 2, na.strings = "999")
  #DON_tab chargeable avec read.table("RGal2012_tout20200401.don",header = TRUE,sep = "\t", skip = 2)
  require("tidyverse")
  require("poppr")
  colnames(DON_tab)[1:2]<-(c("Unit","ind"))
  SSRce_ihc<-VAR_tab%>%
    filter(Unit %in% DON_tab$Unit)
  name_SSRce_ihc<-DON_tab%>%
    filter(Unit %in% SSRce_ihc$Unit)
  SSRce_ihc$Unit<-name_SSRce_ihc$ind
  colnames(SSRce_ihc)[1]<-"ind"
  
  SSRce_ihc_test<-SSRce_ihc[, colSums(SSRce_ihc != 0, na.rm =TRUE) > 0]
  
  colnames(SSRce_ihc_test)<-str_replace(colnames(SSRce_ihc_test),"\\.\\.","\\.999")
  SSRce_gind<-as.genind(tab = as.matrix(SSRce_ihc_test[,-1]),type = "codom", check_ploidy = TRUE)
  
  return(SSRce_gind)
}


allelic_error_computation<-function(genind, ploidy = 2, method = "estimated allele number"){
  # Cette fonction est utilisé pour calculer le taux d'erreur allélique qui est defini tel que
  # Nombre d'erreur / nombre réel d'allèle comparés
  # dependance : poppr, tidyverse
  # Un objet genind doit etre mis en entrée
  locus_error<-as.data.frame(locNames(genind))%>%
    `colnames<-`(.,"Locus")
  locus_error$allelic_error_rate<-0
  locus_error$allelic_mismatch<-0
  locus_error$pair_removed<-0
  locus_error$nAllele_compared<-0
  for (mark in locNames(genind)) {
    print(mark)
    # marker1<-genind$tab %>%
    #   as_tibble()%>%
    #   select(which(str_detect(colnames(genind$tab),mark)))
    
    m<-matrix(data = seq(1:nrow(genind[loc = mark]$tab)),ncol = 2,byrow = TRUE)
    tot_diff<-0
    NA_number<-0
    nAllele<-0
    for (i in seq(1:(nrow(genind[loc = mark]$tab)/2))){
      
      
      if (anyNA(genind[m[i,],loc = mark]$tab)==FALSE) {
        a<-genind[m[i,],loc = mark]%>%
          diss.dist()
        tot_diff<-tot_diff+a
        nAllele<-nAllele + sum(genind[m[i,],loc = mark]$tab)
      }else{NA_number<-NA_number+1}
      
    }
    locus_error$allelic_mismatch[which(locus_error$Locus==mark)]<-tot_diff
    locus_error$pair_removed[which(locus_error$Locus==mark)]<-NA_number
    if (method == "estimated allele number") {
      locus_error$allelic_error_rate[which(locus_error$Locus==mark)]<-tot_diff/((nrow(genind[loc = mark]$tab)-NA_number*2)*ploidy)
    }
    if (method == "sum different allele") {
      locus_error$allelic_error_rate[which(locus_error$Locus==mark)]<-tot_diff/nAllele
    }
    
    locus_error$nAllele_compared[which(locus_error$Locus==mark)]<-nAllele
  }  
  return(locus_error)
}

safe_write_table<-function(tab,file,sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE){
  if (file.exists(paste(file,paste0(str_remove_all(str_sub(Sys.time(),3,10),"-"),".txt"),sep = "_"))==FALSE) {
    write.table(tab,file = paste(file,paste0(str_remove_all(str_sub(Sys.time(),3,10),"-"),".txt"),sep = "_"),col.names = T,row.names = row.names,quote = quote,sep = sep)
    cat("Le data frame", paste(file,paste0(str_remove_all(str_sub(Sys.time(),3,10),"-"),".txt"),sep = "_") ,"a ete enregistre dans le dossier", getwd())
    
  }else{cat("Impossible d'enregistrer car le fichier existe déja")}
}

genotypic_table_to_genalex<-function(geno_tab,file,ploidy = 4){
  # geno_tab : path to the genotypic_table output from FDSTools
  # file : path to where you want to write the new tab
  # Ploidy should be equal to the A parameter use durign the FDSTools calling process
  
  geno_tab$Indiv<-str_remove_all(geno_tab$Indiv,"Polyrosa-")
  geno_tab$pop<-geno_tab$Indiv
  geno_tab$pop[geno_tab$pop %>% str_starts("di-")]<-"Diploids"
  geno_tab$pop[geno_tab$pop %>% str_starts("tetra-")]<-"Tetraploids"  
  geno_tab$pop[geno_tab$pop %>% str_starts("tri-")]<-"Triploids"  
  geno_tab$pop[geno_tab$pop %>% str_starts("penta-")]<-"Pentaploids"
  geno_tab$pop[geno_tab$pop %>% str_starts("hexa-")]<-"Hexaploids"
  geno_tab$pop[geno_tab$pop %>% str_starts("hap-")]<-"Haploid"
  geno_tab$pop[geno_tab$pop %>% str_starts("tneg-")]<-"Tneg"
  geno_tab<-geno_tab%>%mutate_all(~replace(., is.na(.), 0))
  geno_tab[geno_tab == "."]<-"0"
  
  entete<-str_flatten(c(str_flatten(c((ncol(geno_tab)-2)/ploidy,nrow(geno_tab),length(unique(geno_tab$pop))),collapse = ";"),str_flatten(rep(";",ncol(geno_tab)-3))))
  second_line<-str_flatten(rep(";",ncol(geno_tab)-1))
  markers<-colnames(geno_tab)[str_which(colnames(geno_tab),"SSRseq")]
  markers<-str_remove(markers, "SSRseq_Rosa_")
  genalex_colname<-str_flatten(c(c("Ind",";","Pop",";"), str_flatten(markers,collapse = str_flatten(rep(";",ploidy))),str_flatten(rep(";",ploidy-1))))
 
  
  tab<-dplyr::select(geno_tab,c(1,ncol(geno_tab),2:ncol(geno_tab)-1))
  
  write(entete,file,sep = ";",append = FALSE)
  write(second_line,file,sep = ";",append = TRUE,ncolumns = 3)
  write(genalex_colname,file,sep = "",append = TRUE)
  write.table(tab,file,sep = ";",append = TRUE,quote = FALSE,row.names = FALSE,col.names = FALSE)
}

All_number_loc_ind<-function(locus,genind){
  # This function calculates the number of allele per individual per locus
  # locus should be a locus name
  # genind should be a genind object 
  #This function could be run other many loci by using function from the apply family
  all_tab<-as.data.frame(genind@tab)
  n_all<-all_tab%>%
    dplyr::select(contains(locus))%>%
    mutate(sum = rowSums(across(where(is.numeric))))%>%
    dplyr::select(ncol(.))
  colnames(n_all)<-locus
  return(n_all)
}
