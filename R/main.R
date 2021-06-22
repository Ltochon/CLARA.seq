library(TraMineR)
library(cluster)
library(dplyr)
library(doParallel) #librairie pour le paralellisme

#' CLARA clustering
#' @description With the help of TraMineR package, CLARA clustering provide a clustering of big dataset.\cr The main objective is to cluster state sequences with the "LCS" distance calculation method to find the best partition in N clusters.
#'
#' @import TraMineR
#' @import cluster
#' @import dplyr
#' @import doParallel
#'
#' @param data The dataset to use. In case of sequences, use seqdef (from TraMineR package) to create such an object.
#' @param nb_sample The number of subsets to test.
#' @param size_sample The size of each subset
#' @param nb_cluster The number of medoids
#' @param method The calculation method to compute the distance matrix. (See the function seqdist in TraMineR package)
#' @param plot Boolean variable to plot the result of clustering
#' @param find_best_method Method to select the best subset. "Distance" is for the mean distance and "DB" is for Davies-Bouldin value.
#' @param with.diss Boolean if the distance matrix should be returned
#' @param cores Number of cores to use for parallelism
#'
#' @return An object with the data, the medoids id (name of the line), the clustering and the distance matrix
#' @export
#'
#' @examples
#'
#'\dontrun{
#' library(TraMineR)
#' data(biofam)
#' #Basic CLARA Computing
#' my_cluster <- clara_clust(seqdef(biofam), method="LCS")
#'
#' #Improved CLARA Computing
#' my_cluster <- clara_clust(seqdef(biofam), method="LCS", improve_method="Silhouette", nbessai_improved=30)
#' }


clara_clust <- function(data, nb_sample = 100, size_sample = 40 + 2*nb_cluster, nb_cluster = 4, method = "LCS", plot = FALSE, find_best_method = "Distance", with.diss = TRUE, cores = 2){
  message("\nCLARA ALGORITHM Improved\n")
  if(nb_cluster > size_sample){
    stop("Too many cluster requested")
  }
  cl <- cores %>% makeCluster
  registerDoParallel(cl)
  start.time <- proc.time() #debut du processus
  calc_pam <- foreach(loop=1:nb_sample, .packages = c('TraMineR', 'cluster'), .combine = 'c', .multicombine = TRUE, .init = list()) %dopar%{ #on stocke chaque sample
    # for(loop in 1:10){
    data_subset <- data[sample(nrow(data), size_sample),]
    diss <- suppressMessages(seqdist(data_subset, method = method, with.missing = TRUE)) #get matrix dissimilarity
    clustering <- pam(diss,nb_cluster,diss = TRUE, pamonce = 2) #PAM sur la matrice de distance avec fastestPAM
    med <- rownames(data_subset)[clustering$id.med]
    diss <- data.frame()
    for(i in 1:length(med)){
      print("ok")
      if(i == 1){
        diss <- cbind(suppressMessages(seqdist(data, refseq = which(row.names(data)==med[i]), method = method, with.missing = TRUE))) #get matrix dissimilarity
      }
      else{
        diss <- cbind(diss,suppressMessages(seqdist(data, refseq = which(row.names(data)==med[i]), method = method, with.missing = TRUE))) #get matrix dissimilarity

      }
    }
    diss_clustering <- apply(diss, 1,which.min)
    if(find_best_method == "Distance"){
      mean_diss <- mean(apply(diss, 1, min))
    }
    else if(find_best_method == "DB"){
      #####
      seq_obj <- list(seq = data, id.med = med, clusters = diss_clustering, diss = diss)
      list_diam <- vector(mode = "list", length = length(seq_obj$id.med))
      for(i in 1:length(seq_obj$id.med)){ #on stocke chaque sample
        #AMELIORATION CAR DISS RETURNED IN OBJ
        diss <- seq_obj$diss[,i]
        val_diam <- 0
        for(j in 1:length(which(seq_obj$clusters == i))){
          val_diam <- val_diam + (diss[which(seq_obj$clusters == i)[j]])**2
        }
        val_diam <- (val_diam/length(which(seq_obj$clusters == i)))**(1/2)
        list_diam[[i]] <- val_diam
      }
      maximum <- rep(0,length(seq_obj$id.med))
      for(medoid1 in 1:length(seq_obj$id.med)){ #pour chaque sous-groupes
        for(medoid2 in 1:length(seq_obj$id.med)){
          if(medoid1 != medoid2){
            db <- (list_diam[[medoid1]] + list_diam[[medoid2]])/seq_obj$diss[which(rownames(seq_obj$seq) == seq_obj$id.med[medoid1]),medoid2]
            if(db > maximum[medoid1]){
              maximum[medoid1] <- db
            }
          }
        }
      }
      final_db <- mean(maximum)
      mean_diss <- final_db
    }


    #####

    list(mean_diss,diss_clustering,med)
  }
  stopCluster(cl)
  mean_all_diss <- calc_pam[1:length(calc_pam)][seq_along(calc_pam[1:length(calc_pam)]) %% 3 == 1]
  clustering_all_diss <- calc_pam[1:length(calc_pam)][seq_along(calc_pam[1:length(calc_pam)]) %% 3 == 2]
  med_all_diss <- calc_pam[1:length(calc_pam)][seq_along(calc_pam[1:length(calc_pam)]) %% 3 == 0]

  ####diss
  if(with.diss){
    cl <- detectCores() %>% -1 %>% makeCluster
    registerDoParallel(cl)
    calc_diss <- foreach(i=1:length(med_all_diss[[which.min(mean_all_diss)]]), .packages = c('TraMineR', 'cluster'), .combine = 'c', .multicombine = TRUE, .init = list()) %dopar%{ #on stocke chaque sample
      diss <- cbind(suppressMessages(seqdist(data, refseq = which(rownames(data) == med_all_diss[[which.min(mean_all_diss)]][i]), method = method, with.missing = TRUE))) #get matrix dissimilarity
      list(diss)
    }
    stopCluster(cl)
    names(calc_diss) <- 1:length(calc_diss)
    diss <- do.call(cbind,calc_diss)
    bestcluster <- list(seq = data, id.med = med_all_diss[[which.min(mean_all_diss)]], clusters = clustering_all_diss[[which.min(mean_all_diss)]], diss = diss)#création de l'objet à retourner

  }
  else{
    bestcluster <- list(seq = data, id.med = med_all_diss[[which.min(mean_all_diss)]], clusters = clustering_all_diss[[which.min(mean_all_diss)]])#création de l'objet à retourner
  }
  ####diss
  message(paste("\nTable of Sample's Values with", find_best_method, "Method"))
  nb_char = nchar(length(mean_all_diss))#calcul de la longueur du numéro du dernier cluster
  nb_max_med_col1 <- 0
  nb_max_med_col2 <- 0
  for(p in 1:length(med_all_diss)){
    if(p %% 2 == 1){
      if(nb_max_med_col1 < sum(nchar(unlist(med_all_diss[[p]])))){
        nb_max_med_col1 <- sum(nchar(unlist(med_all_diss[[p]])))
      }
    }
    else{
      if(nb_max_med_col2 < sum(nchar(unlist(med_all_diss[[p]])))){
        nb_max_med_col2 <- sum(nchar(unlist(med_all_diss[[p]])))
      }
    }
  }
  nb_precision1  <- 0
  nb_precision2 <- 0
  for(p in 1:length(mean_all_diss)){
    if(p %% 2 == 1){
      if(nb_precision1 < nchar(unlist(mean_all_diss[[p]]))){
        nb_precision1 <- nchar(unlist(mean_all_diss[[p]]))
      }
    }
    else{
      if(nb_precision2 < nchar(unlist(mean_all_diss[[p]]))){
        nb_precision2 <- nchar(unlist(mean_all_diss[[p]]))
      }
    }
  }
  txt = ""
  for(print_dist in 1:length(mean_all_diss)){
    if(print_dist %% 2 == 1){
      #affichage des distances 5 par 5
      txt = paste0(txt,"\nN°",print_dist, strrep(" ",nb_char - nchar(print_dist)), " (Med: ", toString(med_all_diss[[print_dist]]), ")",  strrep(" ",nb_max_med_col1 - sum(nchar(unlist(med_all_diss[[print_dist]])))) ," : ",mean_all_diss[[print_dist]], strrep(" ",nb_precision1 + 2 - nchar(unlist(mean_all_diss[[print_dist]]))))
    }
    else {
      txt = paste0(txt, "|", strrep(" ",2), "N°",print_dist, strrep(" ",nb_char - nchar(print_dist)), " (Med: ", toString(med_all_diss[[print_dist]]), ")",  strrep(" ",nb_max_med_col2 - sum(nchar(unlist(med_all_diss[[print_dist]])))) ," : ",mean_all_diss[[print_dist]], strrep(" ",nb_precision2 + 2 - nchar(unlist(mean_all_diss[[print_dist]]))))
    }

  }
  message(txt)
  #affichage du minimum
  message("\n[>] Minimum Index Value for Sample N°", (which.min(mean_all_diss)),"\n")

  ##########
  #Affichage des graphs
  ##########
  if(plot){
    par(mfrow = c(1,1))
    m <- sapply(seq_along(mean_all_diss),function(x){min(unlist(mean_all_diss)[1:x])})
    index <- c(1,unlist(lapply(1:9,function(x){x*floor(nb_sample/9)})))
    plot(m[index],type = "o", main = paste("Evolution of sample's value with", find_best_method,"method"), xlab = "Iteration Number", ylab = paste(find_best_method ,"value"), col ="blue", pch = 19, lwd = 1)
    seqdplot(data, group = clustering_all_diss[[which.min(mean_all_diss)]], main = "Cluster")
  }
  end.time<-proc.time() #fin du processus
  message("Calculation time : ", (end.time-start.time)[3], " sec.")
  return(bestcluster)
}


#' Davies Bouldin Index
#' @description Implementation of Davies-Bouldin index to evaluate the quality of CLARA Clustering.
#'
#' @import TraMineR
#' @import dplyr
#' @import doParallel
#'
#' @param seq_obj The object generated with CLARA Clustering
#' @param method The calculation method to compute the distance matrix. (See the function seqdist in TraMineR package)
#' @param diss Boolean to express if the parameter diss from CLARA.seq clustering has been returned (Matrix size must be k columns and n rows - see refseq function from TraMineR package)
#' @param plot Boolean variable to plot the result
#' @param cores Number of cores to use for parallelism
#'
#' @return The value of the index
#' @export
#'
#' @examples
#'\dontrun{
#' my_index <- db_index(my_cluster, "LCS", diss = FALSE)
#'}
davies_bouldin <- function(seq_obj, method, diss = TRUE, plot = TRUE, cores = 2){
  message("\nDAVIES-BOULDIN INDEX\n")
  start.time <- proc.time() #debut du processus
  sum_DB <- 0
  res <- 0
  dev_DB = c()
  list_diam <- vector(mode = "list", length = length(seq_obj$id.med))
  if(diss == FALSE){
    cl <- cores %>% makeCluster
    registerDoParallel(cl)
    calc_diss <- foreach(i=1:length(seq_obj$id.med), .packages = c('TraMineR', 'cluster'), .combine = 'c', .multicombine = TRUE, .init = list()) %dopar%{ #on stocke chaque sample
      diss <- cbind(suppressMessages(seqdist(seq_obj$seq, refseq = which(rownames(seq_obj$seq) == seq_obj$id.med[i]), method = method, with.missing = TRUE))) #get matrix dissimilarity
      val_diam <- 0
      for(j in 1:length(which(seq_obj$clusters == i))){
        val_diam <- val_diam + (diss[which(seq_obj$clusters == i)[j]])**2
      }
      val_diam <- (val_diam/length(which(seq_obj$clusters == i)))**(1/2)
      list(val_diam,diss)
    }
    all_diss <- calc_diss[1:length(calc_diss)][seq_along(calc_diss[1:length(calc_diss)]) %% 2 == 0]
    all_diam <- calc_diss[1:length(calc_diss)][seq_along(calc_diss[1:length(calc_diss)]) %% 2 == 1]
    stopCluster(cl)
    maximum <- rep(0,length(seq_obj$id.med))
    for(medoid1 in 1:length(seq_obj$id.med)){ #pour chaque sous-groupes
      for(medoid2 in 1:length(seq_obj$id.med)){
        if(medoid1 != medoid2){
          db <- (all_diam[[medoid1]] + all_diam[[medoid2]])/all_diss[[medoid2]][which(rownames(seq_obj$seq) == seq_obj$id.med[medoid1])]
          if(db > maximum[medoid1]){
            maximum[medoid1] <- db
          }
        }
      }
    }
  }
  else{
    for(i in 1:length(seq_obj$id.med)){ #on stocke chaque sample
      #AMELIORATION CAR DISS RETURNED IN OBJ
      diss <- seq_obj$diss[,i]
      val_diam <- 0
      for(j in 1:length(which(seq_obj$clusters == i))){
        val_diam <- val_diam + (diss[which(seq_obj$clusters == i)[j]])**2
      }
      val_diam <- (val_diam/length(which(seq_obj$clusters == i)))**(1/2)
      list_diam[[i]] <- val_diam
    }
    maximum <- rep(0,length(seq_obj$id.med))
    for(medoid1 in 1:length(seq_obj$id.med)){ #pour chaque sous-groupes
      for(medoid2 in 1:length(seq_obj$id.med)){
        if(medoid1 != medoid2){
          db <- (list_diam[[medoid1]] + list_diam[[medoid2]])/seq_obj$diss[which(rownames(seq_obj$seq) == seq_obj$id.med[medoid1]),medoid2]
          if(db > maximum[medoid1]){
            maximum[medoid1] <- db
          }
        }
      }
    }
  }
  final_db <- mean(maximum)
  db_evolution <- cumsum(maximum) / c(1:length(seq_obj$id.med))
  if(plot){
    plot(db_evolution, main = "Evolution of Davies-Bouldin Index by cluster addition", type = "o", xaxt = "n",xlab = "Cluster number", ylab = "Davies-Bouldin Index", col ="orange", pch = 19, lwd = 1)
    axis(1, at=1:length(seq_obj$id.med), labels= 1:length(seq_obj$id.med))
  }
  message("Value of DB Index for a ", length(seq_obj$id.med), "-clusters : ", final_db)
  end.time<-proc.time() #fin du processus
  message("Calcul time : ", (end.time-start.time)[3], " sec.")
  return(final_db)

}

