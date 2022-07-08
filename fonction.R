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
# print(spe)
dat <- read.table("data_analysis_8pl/GenotypicTable_AlleleDosage_forTest.txt", sep="", header=TRUE)



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
matrix_to_Darwin<-function(genind, file){
  var_RepA<-as.data.frame(tab(genind))%>%
    rownames_to_column()%>%
    rename(Unit= rowname)
  var_RepA[is.na(var_RepA)] = 999
  
  Don_RepA<-as.data.frame(var_RepA$Unit)%>%
    rownames_to_column()
  var_RepA$Unit<-rownames(var_RepA)
  colnames(Don_RepA)<-c("Unit", "Geno")
  line1var<-"@DARwin 5.0 - ALLELIC - 2"
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
  # l'individu retenu sera celui avc le moins de données manquantes
  # un objet genind est retourné
  # Dependences : Polysat, Poppr, tidyverse
  # crée le 06 70 2022 par Clovis
  dist_ind<-diss.dist(genind,percent = TRUE, mat = TRUE)
  # plotdist<-as.data.frame(as.vector(dist_ind))
  # colnames(plotdist)<-"dist"
  # ggplot(data = plotdist)+
  #   geom_density(aes(dist))
  clone_name<-assignClones(dist_ind,threshold = threshold)%>%as.data.frame()%>%
    rownames_to_column()%>%
    reshape::rename(c("rowname"="ind","."="clone" ))%>%
    rownames_to_column()%>%
    mutate(rowname = as.integer(rowname))
  clone_name$missing<-1-propTyped(genind,by="ind") 

  ind_to_kept<-clone_name%>%
    group_by(clone)%>%
    slice_min(missing,n = 1,with_ties = FALSE)
  
  
  mdata_cl_rmvd<- genind[i = ind_to_kept$rowname, drop = TRUE]
  return(mdata_cl_rmvd)
}


recall_genind <- function(allele_tab,correspondance, prs_abs = FALSE) {
  a<-genind2df(allele_tab,sep = "_")
  i<-2
  # allele_tab : un objet genind
  # Un data.frame des correspondances entre code allélique
  # Cette fonction permet de passer de changer le code allélique des individus 
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
    select(-pop)
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
