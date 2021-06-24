#' CLARA clustering
#' @description With the help of TraMineR package, CLARA clustering provide a clustering of big dataset.\cr The main objective is to cluster state sequences with the "LCS" distance calculation method to find the best partition in N clusters.
#'
#' @import TraMineR
#' @import cluster
#' @import dplyr
#' @import doParallel
#' @import parallel
#' @import foreach
#'
#' @param data The dataset to use. In case of sequences, use seqdef (from TraMineR package) to create such an object.
#' @param nb_sample The number of subsets to test.
#' @param size_sample The size of each subset
#' @param nb_cluster The number of medoids
#' @param distargs List with method parameters to apply. (See the function seqdist in TraMineR package)
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
#' #creating sequences
#' library(TraMineR)
#' data(mvad)
#' mvad.labels <- c("employment", "further education", "higher education","joblessness", "school", "training")
#' mvad.scode <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' mvad.seq <- seqdef(mvad, 17:86, states = mvad.scode,abels = mvad.labels, xtstep = 6)
#'
#' #CLARA Clustering
#' my_cluster <- clara_clust(mvad.seq,nb_cluster = 4, nb_sample = 10, size_sample = 20, with.diss = TRUE)
#'
#' #CLARA Clustering with Davies-Bouldin Method
#' my_cluster <- clara_clust(mvad.seq,nb_cluster = 4, nb_sample = 10, size_sample = 20, with.diss = TRUE, find_best_method = "DB")


clara_clust <- function(data, nb_sample = 100, size_sample = 40 + 2*nb_cluster, nb_cluster = 4, distargs = list(method = "LCS"), plot = FALSE, find_best_method = "Distance", with.diss = TRUE, cores = detectCores()-1){
  message("\nCLARA ALGORITHM Improved\n")
  if(nb_cluster > size_sample){
    stop("Too many cluster requested")
  }
  cl <- cores %>% makeCluster
  registerDoParallel(cl)
  start.time <- proc.time() #debut du processus
  calc_pam <- foreach(loop=1:nb_sample, .packages = c('TraMineR', 'cluster'), .combine = 'c', .multicombine = TRUE, .init = list()) %dopar%{ #on stocke chaque sample
    # for(loop in 1:10){
    distargs$refseq <- NULL
    data_subset <- data[sample(nrow(data), size_sample),]
    distargs$seqdata <- data_subset
    suppressMessages(diss <- do.call(seqdist, distargs))
    clustering <- pam(diss,nb_cluster,diss = TRUE, pamonce = 2) #PAM sur la matrice de distance avec fastestPAM
    med <- rownames(data_subset)[clustering$id.med]
    diss <- data.frame()
    distargs$seqdata <- data
    for(i in 1:length(med)){
      distargs$refseq <- which(row.names(data)==med[i])
      if(i == 1){
        # diss <- cbind(suppressMessages(seqdist(data, refseq = which(row.names(data)==med[i]), method = method, with.missing = TRUE))) #get matrix dissimilarity
        diss <- cbind(suppressMessages(diss <- do.call(seqdist, distargs)))
      }
      else{
        # diss <- cbind(diss,suppressMessages(seqdist(data, refseq = which(row.names(data)==med[i]), method = method, with.missing = TRUE))) #get matrix dissimilarity
        diss <- cbind(diss,suppressMessages(diss <- do.call(seqdist, distargs)))
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
    cl <- cores %>% makeCluster
    registerDoParallel(cl)
    calc_diss <- foreach(i=1:length(med_all_diss[[which.min(mean_all_diss)]]), .packages = c('TraMineR', 'cluster'), .combine = 'c', .multicombine = TRUE, .init = list()) %dopar%{ #on stocke chaque sample
      # diss <- cbind(suppressMessages(seqdist(data, refseq = which(rownames(data) == med_all_diss[[which.min(mean_all_diss)]][i]), method = method, with.missing = TRUE))) #get matrix dissimilarity
      distargs$seqdata <- data
      distargs$refseq <- which(rownames(data) == med_all_diss[[which.min(mean_all_diss)]][i])
      diss <- cbind(suppressMessages(diss <- do.call(seqdist, distargs)))
      list(diss)
    }
    stopCluster(cl)
    names(calc_diss) <- 1:length(calc_diss)
    diss <- do.call(cbind,calc_diss)
    bestcluster <- list(seq = data, id.med = med_all_diss[[which.min(mean_all_diss)]], clusters = clustering_all_diss[[which.min(mean_all_diss)]], diss = diss, evol.diss = sapply(seq_along(mean_all_diss),function(x){min(unlist(mean_all_diss)[1:x])}), plot = function(){plot(bestcluster$evol.diss[c(1,unlist(lapply(1:9,function(x){x*floor(nb_sample/9)})))],type = "o", main = paste("Evolution of sample's value with", find_best_method,"method"), xlab = "Iteration Number", ylab = paste(find_best_method ,"value"), col ="blue", pch = 19, lwd = 1)})#création de l'objet à retourner

  }
  else{
    bestcluster <- list(seq = data, id.med = med_all_diss[[which.min(mean_all_diss)]], clusters = clustering_all_diss[[which.min(mean_all_diss)]], evol.diss = sapply(seq_along(mean_all_diss),function(x){min(unlist(mean_all_diss)[1:x])}), plot = function(){plot(bestcluster$evol.diss[c(1,unlist(lapply(1:9,function(x){x*floor(nb_sample/9)})))],type = "o", main = paste("Evolution of sample's value with", find_best_method,"method"), xlab = "Iteration Number", ylab = paste(find_best_method ,"value"), col ="blue", pch = 19, lwd = 1)})#création de l'objet à retourner
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
    plot(m[index],type = "o", main = paste("Evolution of sample's value with", find_best_method,"method"), xlab = "Iteration Number", ylab = paste(find_best_method ,"value"), col ="blue", pch = 19, lwd = 1, xaxt = "n")
    axis(1, at = 1:10, labels = index)
    # seqdplot(data, group = clustering_all_diss[[which.min(mean_all_diss)]], main = "Cluster")
  }
  end.time<-proc.time() #fin du processus
  message("Calculation time : ", (end.time-start.time)[3], " sec.")
  class(bestcluster) <- c("clara", "partition")
  return(bestcluster)
}


#' CLARANS clustering
#' @description With the help of TraMineR package, CLARANS clustering provide a clustering of big dataset.\cr The main objective is to cluster state sequences with the "LCS" distance calculation method to find the best partition in N clusters. \cr \cr WARNING : this function is less efficient than cLARA.
#'
#' @import TraMineR
#' @import cluster
#' @import dplyr
#' @import doParallel
#'
#' @param data The dataset to use. In case of sequences, use seqdef (from TraMineR package) to create such an object.
#' @param nb_cluster The number of medoids
#' @param distargs List with method parameters to apply. (See the function seqdist in TraMineR package)
#' @param maxneighbours Number of neighbours to explore to find a better clustering
#' @param numlocal Number of initialisation of the starting medoids
#' @param plot Boolean variable to plot the research convergence
#' @param cores Number of cores to use for parallelism
#'
#' @return An object with the data, the medoids id (name of the line), the clustering and the distance matrix
#' @export
#'
#' @examples
#'
#' #creating sequences
#' library(TraMineR)
#' data(mvad)
#' mvad.labels <- c("employment", "further education", "higher education","joblessness", "school", "training")
#' mvad.scode <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' mvad.seq <- seqdef(mvad, 17:86, states = mvad.scode,abels = mvad.labels, xtstep = 6)
#'
#' #CLARANS Clustering
#' my_cluster <- clarans_clust(mvad.seq,nb_cluster = 4, maxneighbours = 20, numlocal = 4, plot = TRUE)

clarans_clust <- function(data, nb_cluster, distargs = list(method = "LCS"), maxneighbours, numlocal, plot = FALSE, cores = detectCores()-1){
  message("\nCLARANS ALGORITHM\n")
  start.time<-proc.time()
  '%ni%' = Negate('%in%')
  mincost <- Inf
  cl <- cores %>% makeCluster
  registerDoParallel(cl)
  final_list = list()
  for(i in 1:numlocal){
    medoids <- sample(1:length(data[,1]),nb_cluster) #choose random med
    while(anyDuplicated(data[medoids,]) > 0){ #check if 2 meds have same value
      medoids <- sample(1:length(data[,1]),nb_cluster)
    }
    calc_diss <- foreach(j=1:nb_cluster, .packages = c('TraMineR', 'cluster'), .combine = 'c', .multicombine = TRUE, .init = list()) %dopar%{
      distargs$refseq <- medoids[j]
      distargs$seqdata <- data
      # diss_current <- cbind(suppressMessages(seqdist(data, refseq = medoids[j], method = method, with.missing = TRUE))) #get matrix dissimilarity
      cbind(suppressMessages(diss_current <- do.call(seqdist, distargs)))
      list(diss_current)
    }
    diss_current <- as.data.frame(do.call(cbind, calc_diss))
    diss_test <- diss_current
    current_sum <- sum(apply(diss_current, 1,min))
    j = 1
    calc_bestseq <- foreach(j = 1:maxneighbours, .packages = c('TraMineR', 'cluster'), .combine = 'c', .multicombine = TRUE, .init = list(list(),list(),list())) %dopar%{
      next_iter = 0
      while(next_iter == 0){
        out_medoid <- sample(1:nb_cluster,1)
        candidate <- sample(setdiff(1:length(data[,1]),medoids),1)
        while(nrow(merge(data[candidate,],data[medoids,]))>1){ #test if candidate has same value than a medoid
          candidate <- sample(setdiff(1:length(data[,1]),medoids),1)
        }
        distargs$seqdata <- data
        distargs$refseq <- candidate
        suppressMessages(diss_candidate <- do.call(seqdist, distargs))
        # diss_candidate <- suppressMessages(seqdist(data, refseq = candidate, method = method, with.missing = TRUE)) #get matrix dissimilarity
        diss_test[,out_medoid] <- diss_candidate
        test_sum <- sum(apply(diss_test, 1,min))
        if(test_sum < current_sum){
          medoids[out_medoid] <- candidate
          diss_current <- diss_test
          current_sum <- test_sum
        }
        else{
          next_iter = 1
        }
      }
      if(j == maxneighbours){
        list(current_sum,medoids,diss_current)
      }
    }
    final_list <- c(final_list,calc_bestseq[4:length(calc_bestseq)])
  }
  stopCluster(cl)
  all_diss <- final_list[seq_along(final_list) %% 3 == 0]
  all_medoids <- final_list[seq_along(final_list) %% 3 == 2]
  all_sum <- final_list[seq_along(final_list) %% 3 == 1]
  #création de la classe
  bestcluster <- list(seq = data, id.med = all_medoids[[which.min(all_sum)]], clusters = apply(all_diss[[which.min(all_sum)]], 1,which.min), diss = all_diss[[which.min(all_sum)]])#création de l'objet à retourner

  message("Table of Iteration's Distance")
  nb_char <- nchar(length(all_sum))#calcul de la longueur du numéro du dernier cluster
  nb_precision <- 17 #on définit le nombre de caractères pour une distance avec la précision maximale (14 décimales)
  nb_max_med <- nchar(length(data[,1])) * nb_cluster + (nb_cluster-1)*2
  txt = ""
  t <- c()
  for(print_dist in 1:length(all_sum)){
    t <- c(t,mean(apply(all_diss[[print_dist]], 1, min)))
    #affichage des distances 5 par 5
    if(print_dist %% 2 == 1){
      txt = paste0(txt,"\nN°",print_dist, strrep(" ",nb_char - nchar(print_dist)), " (Med: ", toString(all_medoids[[print_dist]]), ")", strrep(" ",nb_max_med - nchar(toString(all_medoids[[print_dist]]))) ," : ",mean(apply(all_diss[[print_dist]], 1, min)), strrep(" ",nb_precision - nchar(mean(apply(all_diss[[print_dist]], 1, min)))))
    }
    else {
      txt = paste0(txt,"  ", "N°",print_dist, strrep(" ",nb_char - nchar(print_dist)),  " (Med: ", toString(all_medoids[[print_dist]]), ")", strrep(" ",nb_max_med - nchar(toString(all_medoids[[print_dist]]))) , " : ",mean(apply(all_diss[[print_dist]], 1, min)), strrep(" ",nb_precision - nchar(mean(apply(all_diss[[print_dist]], 1, min)))))
    }
  }
  message(txt)
  #affichage du minimum
  message("\n[>] Minimum Distance for iteration N°", (which.min(all_sum)),"\n")

  ##########
  #Affichage des graphs
  ##########
  if(plot){
    par(mfrow = c(1,1))
    #affichage du graph de variation des distances pour un nombre condensé de sample
    plot(sapply(seq_along(all_diss), function(x){min(t[1:x])}), type = "o", main = "Evolution of the mean distance", xlab = "NumLocal number", ylab = "Mean distance value", col ="blue", pch = 19, lwd = 1)
    # seqdplot(data, apply(all_diss[[which.min(all_sum)]], 1,which.min), main = "Cluster")
  }
  end.time<-proc.time() #fin du processus
  message("Calculation time : ", (end.time-start.time)[3], " sec.")
  class(bestcluster) <- c("clarans", "partition")
  return(bestcluster)
}

#' FUZZY-CLARA clustering
#' @description With the help of TraMineR package, FUZZY-CLARA clustering provide a clustering of big dataset.\cr The main objective is to cluster state sequences with the "LCS" distance calculation method to find the best partition in N clusters. \cr This function is a mix between CLARA and the ROBUST FUZZY C-MEDOIDS. \cr \cr WARNING : This function is not finished yet !
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
#' @param distargs List with method parameters to apply. (See the function seqdist in TraMineR package)
#' @param fuzzyfier Value of the fuzzifier (default is 2, which is the traditionnal value)
#' @param p Number of candidate to test to be a better medoid
#' @param threshold Variable to exclude outliers, whose values are greater than threshold
#' @param max_iter Number of maximal iteration to do to find the set of medoids
#' @param noise Small value to avoid divisions by 0 error
#' @param plot Boolean variable to plot the research convergence
#' @param cores Number of cores to use for parallelism
#'
#' @return An object with the data, the medoids id (name of the line), the clustering and the distance matrix
#' @export
#'
#' @examples
#'
#' #creating sequences
#' library(TraMineR)
#' data(mvad)
#' mvad.labels <- c("employment", "further education", "higher education","joblessness", "school", "training")
#' mvad.scode <- c("EM", "FE", "HE", "JL", "SC", "TR")
#' mvad.seq <- seqdef(mvad, 17:86, states = mvad.scode,abels = mvad.labels, xtstep = 6)
#'
#' #CLARA-FUzZY Clustering
#' my_cluster <- fuzzy_clust(mvad.seq,nb_sample = 14, size_sample = 50, plot = TRUE, threshold = 7, max_iter = 10, p=5)
#'
fuzzy_clust <- function(data, nb_sample = 100, size_sample = 40 + 2*nb_cluster, nb_cluster = 4, distargs = list(method = "LCS"), fuzzyfier = 2, p = 5, threshold = 10, max_iter = 10, noise = 0.5, plot = FALSE, cores = detectCores()-1){
  if(nb_cluster > size_sample){
    stop("Too many cluster requested")
  }
  message("ROBUST FUZZY C-MEDOIDS ALGORITHM")
  cl <- cores %>% makeCluster
  registerDoParallel(cl)
  start.time <- proc.time() #debut du processus
  calc_fuzzy <- foreach(loop=1:nb_sample, .packages = c('TraMineR'), .combine = 'c', .multicombine = TRUE, .init = list()) %dopar%{ #on stocke chaque sample
    # for(loop in nb_sample){
    data_subset <- data[sample(nrow(data), size_sample),]
    distargs$seqdata <- data_subset
    suppressMessages(diss <- do.call(seqdist, distargs))
    # diss <- suppressMessages(seqdist(data_subset, method = method, with.missing = TRUE)) #get matrix dissimilarity
    med = rownames(data_subset)[sample(size_sample, nb_cluster)]
    unique <- 0
    while(unique < 1){
      unique <- 1
      for(m in 1:nb_cluster){
        for(n in 1:nb_cluster){
          if(m != n){
            if(seqcomp(data_subset[match(med[m],rownames(data_subset)),],data_subset[match(med[n],rownames(data_subset)),])){
              tab <- match(med,rownames(data_subset))
              med[m] <- rownames(data_subset)[sample(c(1:length(data_subset[,1]))[-c(match(med,rownames(data_subset)))],1)]
              unique <- 0
            }
          }
        }
      }
    }
    med_old = c(0)
    med_old2 = c(0)
    iter = 0
    stop = 0
    while(stop == 0){
      not_valid <- c()
      tab_harm <- c()
      membership <- data.frame()
      for(i in 1:length(data_subset[,1])){
        u_tot <- c()
        d <- diss[match(med,rownames(data_subset)),i]
        for(j in 1:nb_cluster){
          u_up <- (1/(d[j]+noise))^(1/(fuzzyfier-1))
          u_down <- sum((1/(d+noise))^(1/(fuzzyfier-1)))
          u_tot <- c(u_tot,u_up/u_down)
        }
        tab_harm <- c(tab_harm,(sum((1/(d+noise))^(1/(fuzzyfier-1))))^(1-fuzzyfier))
        if((sum((1/(d+noise))^(1/(fuzzyfier-1))))^(1-fuzzyfier) < threshold){
          membership <- rbind(membership, u_tot)
        }
        else{
          not_valid <- c(not_valid,i)
        }
      }
      if(is.null(not_valid)){
        rownames(membership) <- make.names(rownames(data_subset))
      }
      else{
        rownames(membership) <- make.names(rownames(data_subset)[-not_valid])
      }
      med_old2 <- med_old
      med_old <- med
      for(i in 1:nb_cluster){
        xpi <- order(membership[,i], decreasing = TRUE)[1:(p+nb_cluster)]
        if(sum(match(med,substring(row.names(membership),2)) %in% xpi) > 0){
          index_na <- match(match(med,substring(row.names(membership),2)),xpi)
          index <- index_na[!is.na(index_na)]
          nb_index_delete <- nb_cluster - sum(match(med,substring(row.names(membership),2)) %in% xpi)
          xpi <- xpi[-index]
          xpi <- xpi[1:(length(xpi)-nb_index_delete)]
        }
        else{
          xpi <- xpi[1:(length(xpi)-nb_cluster)]
        }
        xpi <- xpi[!is.na(xpi)]
        q_min = Inf
        q_index_min = med[i]
        if(length(xpi) != 0){
          for(j in 1:length(xpi)){
            a <- data_subset[which(rownames(data_subset) == substring(rownames(membership)[xpi[j]],2)),]
            b_index <- match(med,rownames(data_subset))
            double <- 0
            for(t in 1:length(b_index)){
              if(seqcomp(a,data_subset[b_index[t],])){
                double <- 1
              }
            }
            if(double < 1){
              q = (membership[xpi[j],i]^fuzzyfier)*sum(diss[match(substring(row.names(membership[xpi[j],]),2),rownames(data_subset)),])
              if(q < q_min){
                q_min <- q
                q_index_min = substring(rownames(membership)[xpi[j]],2)
              }
            }
          }
          med[i] = q_index_min
        }
      }
      iter = iter + 1
      if(setequal(med_old2,med) || iter >= max_iter){
        stop = 1
      }
    }
    ##############################################
    #############################################
    diss_final <- data.frame()
    distargs$seqdata <- data
    for(i in 1:length(med)){
      distargs$refseq <- which(row.names(data)==med[i])
      if(i == 1){
        diss_final <- cbind(suppressMessages(diss_final <- do.call(seqdist, distargs)))
        # diss_final <- cbind(suppressMessages(seqdist(data, refseq = which(row.names(data)==med[i]), method = method, with.missing = TRUE))) #get matrix dissimilarity
      }
      else{
        diss_final <- cbind(diss_final,suppressMessages(diss_final <- do.call(seqdist, distargs)))
        # diss_final <- cbind(diss_final,suppressMessages(seqdist(data, refseq = which(row.names(data)==med[i]), method = method, with.missing = TRUE))) #get matrix dissimilarity

      }
    }
    membership <- data.frame()
    for(i in 1:length(data[,1])){
      u_tot <- c()
      d <- diss_final[i,]
      for(j in 1:nb_cluster){
        u_up <- (1/(d[j]+noise))^(1/(fuzzyfier-1))
        u_down <- sum((1/(d+noise))^(1/(fuzzyfier-1)))
        u_tot <- c(u_tot,u_up/u_down)
      }
      membership <- rbind(membership, u_tot)
    }
    q <- 0
    for(i in 1:nb_cluster){
      q <- q + (membership[which(rownames(data) == med[i]),i]^fuzzyfier)*sum(diss_final[i,])
    }
    q <- q / nb_cluster
    mean_diss <-mean(apply(membership, 1,max))
    diss_clustering <- apply(membership, 1,which.max)
    list(diss_final,mean_diss,diss_clustering,med,tab_harm,iter,membership,q)
  }
  stopCluster(cl)
  all_diss <- calc_fuzzy[1:length(calc_fuzzy)][seq_along(calc_fuzzy[1:length(calc_fuzzy)]) %% 8 == 1]
  mean_all_diss <- calc_fuzzy[1:length(calc_fuzzy)][seq_along(calc_fuzzy[1:length(calc_fuzzy)]) %% 8 == 2]
  clustering_all_diss <- calc_fuzzy[1:length(calc_fuzzy)][seq_along(calc_fuzzy[1:length(calc_fuzzy)]) %% 8 == 3]
  med_all_diss <- calc_fuzzy[1:length(calc_fuzzy)][seq_along(calc_fuzzy[1:length(calc_fuzzy)]) %% 8 == 4]
  harm_all <- calc_fuzzy[1:length(calc_fuzzy)][seq_along(calc_fuzzy[1:length(calc_fuzzy)]) %% 8 == 5]
  iter_all <- calc_fuzzy[1:length(calc_fuzzy)][seq_along(calc_fuzzy[1:length(calc_fuzzy)]) %% 8 == 6]
  all_memb <- calc_fuzzy[1:length(calc_fuzzy)][seq_along(calc_fuzzy[1:length(calc_fuzzy)]) %% 8 == 7]
  all_q <- calc_fuzzy[1:length(calc_fuzzy)][seq_along(calc_fuzzy[1:length(calc_fuzzy)]) %% 8 == 0]
  message("\nTable of objective function values")
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
      if(nb_precision1 < nchar(unlist(all_q[[p]]))){
        nb_precision1 <- nchar(unlist(all_q[[p]]))
      }
    }
    else{
      if(nb_precision2 < nchar(unlist(all_q[[p]]))){
        nb_precision2 <- nchar(unlist(all_q[[p]]))
      }
    }
  }
  txt = ""
  for(print_dist in 1:length(mean_all_diss)){
    if(print_dist %% 2 == 1){
      #affichage des distances 5 par 5
      txt = paste0(txt,"\nN°",print_dist, strrep(" ",nb_char - nchar(print_dist)), " (Med: ", toString(med_all_diss[[print_dist]]), ")",  strrep(" ",nb_max_med_col1 - sum(nchar(unlist(med_all_diss[[print_dist]])))) ," : ",all_q[[print_dist]], strrep(" ",nb_precision1 + 2 - nchar(unlist(all_q[[print_dist]]))))
    }
    else {
      txt = paste0(txt, "|", strrep(" ",2), "N°",print_dist, strrep(" ",nb_char - nchar(print_dist)), " (Med: ", toString(med_all_diss[[print_dist]]), ")",  strrep(" ",nb_max_med_col2 - sum(nchar(unlist(med_all_diss[[print_dist]])))) ," : ",all_q[[print_dist]], strrep(" ",nb_precision2 + 2 - nchar(unlist(all_q[[print_dist]]))))
    }

  }
  message(txt)
  #affichage du minimum
  message("\n[>] Minimum Index for Sample N°", (which.min(all_q)),"\n")
  #création de la classe
  bestcluster <- list(seq = data, id.med = med_all_diss[[which.min(all_q)]], clusters = clustering_all_diss[[which.min(all_q)]], diss = all_diss[[which.min(all_q)]], harm_value = harm_all[[which.min(all_q)]], nb_iter = iter_all[[which.min(all_q)]], membership = all_memb[[which.min(all_q)]], q = all_q)#création de l'objet à retourner

  ##########
  #Affichage des graphs
  ##########
  if(plot){
    par(mfrow = c(1,1))
    #affichage du graph de variation des distances pour un nombre condensé de sample
    plot(sapply(seq_along(all_q), function(x){min(unlist(all_q)[1:x])}), type = "o", main = "Evolution of the objective function", xlab = "Iteration number", ylab = "Index value", col ="blue", pch = 19, lwd = 1)
    # seqdplot(data, clustering_all_diss[[which.min(all_q)]], main = "Cluster")
  }
  end.time<-proc.time() #fin du processus
  message("Calcul time : ", (end.time-start.time)[3], " sec.")
  class(bestcluster) <- c("clara-fuzzy", "partition")
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
#' @param distargs List with method parameters to apply. (See the function seqdist in TraMineR package)
#' @param diss Boolean to express if the parameter diss from CLARA.seq clustering has been returned (Matrix size must be k columns and n rows - see refseq function from TraMineR package)
#' @param plot Boolean variable to plot the research convergence
#' @param cores Number of cores to use for parallelism
#'
#' @return The value of the index
#' @export
#'
#' @examples
#'\dontrun{
#' my_index <- davies_bouldin(my_cluster)
#'}

davies_bouldin <- function(seq_obj, distargs = list(method = "LCS"), diss = TRUE, plot = TRUE, cores = detectCores()-1){
  message(paste("\nDAVIES-BOULDIN INDEX for", class(seq_obj)[1],"Clustering\n"))
  if(!(is(seq_obj,c("clara", "partition")) || is(seq_obj,c("clarans", "partition")) || is(seq_obj,c("clara-fuzzy", "partition")))){
    stop("Seq_obj must come from clarans_clust, clara_clust or fuzzy_clust")
  }
  start.time <- proc.time() #debut du processus
  sum_DB <- 0
  res <- 0
  dev_DB = c()
  list_diam <- vector(mode = "list", length = length(seq_obj$id.med))
  if(diss == FALSE){
    cl <- cores %>% makeCluster
    registerDoParallel(cl)
    suppressWarnings(distargs$seqdata <- seq_obj$seq)
    calc_diss <- foreach(i=1:length(seq_obj$id.med), .packages = c('TraMineR', 'cluster'), .combine = 'c', .multicombine = TRUE, .init = list()) %dopar%{ #on stocke chaque sample
      distargs$refseq <- which(rownames(seq_obj$seq) == seq_obj$id.med[i])
      diss <- suppressMessages(diss <- do.call(seqdist, distargs))
      # diss <- cbind(suppressMessages(seqdist(seq_obj$seq, refseq = which(rownames(seq_obj$seq) == seq_obj$id.med[i]), method = method, with.missing = TRUE))) #get matrix dissimilarity
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

