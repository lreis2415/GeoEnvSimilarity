#similarity-based inference


#' @description Compute similarity of two locations at the covariate level
#' The equation used in this function is explained in the work of 
#' Zhu et al.(2015, https://doi.org/10.1111/ejss.12244)
#'
#' @param v1 the covariate value of the sample
#' @param v2 the covariate value of the unpredicted location
#' @param type the data type of the covariate, either "CONTINOUS" or "CATEGORICAL"
#' @param sd the standard deviation of the covariate in the study area
#' @param mean the mean of the covariate in the study area
#' @return the similarity of two locations at the covariate level
#' @example 
#' single_simi(3,4,"CONTINUOUS",1,4)
single_simi<-function(v1,v2,type,sd,mean){
  if(type=="CONTINUOUS"){
    v=exp(-(v1-v2)^2*0.5/sd^4*(v1^2+sd^2+mean^2-2*v1*mean)) #ipsm, Zhu et al., 2015
    #v = exp(((v1-v2)/sd)^2*log(0.5)) #range rule
    
    if(v>1) {
      tmp=1 
    } else if(v<0) {
      tmp=0
    } else {
      tmp=v
    }
    return(tmp)
  }
  else {
    if(smpEnv[k]==locEnv[k]) { return(1) }
    else { return(0) }
  }
}

#' @description Compute similarity of two locations
#' The equation used in this function is explained in the work of 
#' Zhu et al.(2015, https://doi.org/10.1111/ejss.12244)
#'
#' @param smpEnv a vector of covariate values of the sample
#' @param locEnv a vector of covariate values of the unpredicted location
#' @param types a vector of the covariates' data types, either "CONTINOUS" or "CATEGORICAL"
#' @param sds a vector of standard deviations of the covariates in the study area
#' @param means a vector of the mean values of the covariates in the study area
#' @param integrate.method the method used for the integration of similarity at covariate level to location level, can be 'MIN', 'MEAN' or 'WEIGHT'(weighted mean), default is 'MIN'
#' @param weight if integrate.method is 'WEIGHT', a vector of weights of different covatiates need to be provided
#' @return the integrated similarity of two locations
loc_simi<-function(smpEnv,locEnv,types,sds,means,integrate.method='MIN', weight){
  cov_num = min(length(smpEnv),length(locEnv),
                length(types),length(sds),length(means))
  simi_tmp=c()
  for(i in 1:cov_num){
    simi_tmp[i]=single_simi(as.numeric(smpEnv[i]),as.numeric(locEnv[i]),
                             types[i],sds[i],means[i])
     
  }
  if(integrate.method == "MEAN"){
    return(mean(simi_tmp))
  } else if(integrate.method == "WEIGHT"){
    if(length(weight)>=cov_num){
      return((simi_tmp*weight)/sum(weight))
    } else {
      stop("The number of weights is inconsistent with the number of covariates.")
    }
  } else {
    return(min(simi_tmp))
  }
}


#' @description Inference based on the similarity between samples and unpredicted locations
#' The method is explained in the work of Zhu et al.(2015,2018 https://doi.org/10.1111/ejss.12244 https://doi.org/10.1080/19475683.2018.1534890)
#'
#' @param samples_cov a numeric list of covariate values of all samples, with column being different covariates and row being different samples
#' @param samples_val a numeric vector of the property values at samples
#' @param locations_cov the covariate values of all unpredicted locations, , with column being different covariates and row being different samples
#' @param types a numeric vector of the covariates' data types, either "CONTINOUS" or "CATEGORICAL"
#' @param threshold a value between 0 and 1. For each unpredicted location, only samples with similarity measures higher than the threshold will be used for prediction. default is 0.
#' @param integrate.method the method used for the integration of similarity at covariate level to location level, can be 'MIN', 'MEAN' or 'WEIGHT'(weighted mean), default is 'MIN'
#' @param weight if integrate.method is 'WEIGHT', a vector of weights of different covatiates need to be provided
#' @return uncertainty: the vector of uncertainty measures at all unpredicted locations
#' @return predicted: the vector of predicted values at all unpredicted locations0
solim<-function(samples_cov,samples_val,locations_cov,types,threshold = 0,
               integrate.method='MIN', weight){
  uncertainty = c()
  predicted = c()
  if(ncol(samples_cov)!=ncol(locations_cov)){
    warning("The covariate number of samples and that of unpredicted locations are different.")
  }
  if(length(samples_val)!=nrow(samples_cov)){
    stop("The number of samples indicated by samples_val and samples_cov are inconsistent.")
  }
  sds<-c()
  means<-c()
  for(k in 1:min(ncol(samples_cov),ncol(locations_cov))){
    sds[k] = sd(append(samples_cov[,k],locations_cov[,k]))
    means[k] = mean(append(samples_cov[,k],locations_cov[,k]))
  }
  for(i in 1:nrow(locations_cov)){
    max_simi=0
    simi_vals=c()
    for(j in 1:nrow(samples_cov)){
      simi_vals[j] = loc_simi(samples_cov[j,],locations_cov[i,],
                              types,sds,means,integrate.method,weight)
    }
    uncertainty[i] = 1-max(simi_vals)
    predicted[i] = sum(samples_val[simi_vals>threshold]*simi_vals[simi_vals>threshold])/sum(simi_vals[simi_vals>threshold])
  }
  result = c()
  result$uncertainty = uncertainty
  result$predicted = predicted
  return(result)
}

#' @description Inference based on the saptial distance and similarity between samples and unpredicted locations
#' The method is explained in the work of Qin et al.(2021, https://doi.org/10.1016/S1002-0160(20)60016-9)
#'
#' @param samples_cov a numeric list of covariate values of all samples, with column being different covariates and row being different samples
#' @param samples_coord a numeric list of coordinates (x,y) of all samples, with two column being x and y, respectively, and row being different samples
#' @param samples_val a numeric vector of the property values at samples
#' @param locations_cov the covariate values of all unpredicted locations, with column being different covariates and row being different unpredicted locations
#' @param samples_coord a numeric list of coordinates (x,y) of all unpredicted locations, with two column being x and y, respectively, and row being different unpredicted locations
#' @param types a vector of the covariates' data types, the element in the vector should be either "CONTINOUS" or "CATEGORICAL"
#' @param threshold a value between 0 and 1. For each unpredicted location, only samples with similarity measures higher than the threshold will be used for prediction. default is 0.
#' @param r a positive real number, default 0.5; the power parameter for calculating the weight of spatial distance. SoLIM-IDW(r=0): the original SoLIM
#' @param integrate.method 'MIN', 'MEAN' or 'WEIGHT'(weighted mean), default is 'MIN'
#' @param weight if integrate.method is 'WEIGHT', a vector of weights of different covatiates need to be provided
#' @return uncertainty: the vector of uncertainty measures at all unpredicted locations
#' @return predicted: the vector of predicted values at all unpredicted locations0
solim.idw<-function(samples_cov,samples_coord,samples_val,locations_cov,locations_coord,
                    types,threshold,r=0.5,integrate.method='MIN', weight){
  uncertainty = c()
  predicted = c()
  if(ncol(samples_cov)!=ncol(locations_cov)){
    warning("The covariate number of samples and that of unpredicted locations are different.")
  }
  if(length(samples_val)!=nrow(samples_cov)){
    stop("The number of samples indicated by samples_val and samples_cov are inconsistent.")
  }
  if(ncol(samples_coord)!=2 | nrow(samples_coord)!=nrow(samples_cov)){
    stop("The parameter samples_coords (coordinates of samples) are not correct")
  }
  if(ncol(locations_coord)!=2 | nrow(locations_coord)!=nrow(locations_cov)){
    stop("The parameter locations_coord (coordinates of unpredicted locations) are not correct")
  }
  sds<-c()
  means<-c()
  for(k in 1:min(ncol(samples_cov),ncol(locations_cov))){
    sds[k] = sd(append(samples_cov[,k],locations_cov[,k]))
    means[k] = mean(append(samples_cov[,k],locations_cov[,k]))
  }
  for(i in 1:nrow(locations_cov)){
    max_simi=0
    simi_vals=c()
    dist_weights = c()
    for(j in 1:nrow(samples_cov)){
      simi_vals[j] = loc_simi(samples_cov[j,],locations_cov[i,],
                              types,sds,means,integrate.method,weight)
      dist_weights[j] = 1/(sqrt((samples_coord[j,1]-locations_cov[i,1])^2+
                        (samples_coord[j,2]-locations_cov[i,2])^2))^r
    }
    uncertainty[i] = 1-max(simi_vals)
    predicted[i] = sum(samples_val[simi_vals>threshold]*simi_vals[simi_vals>threshold]
                       *dist_weights[simi_vals>threshold])/sum(simi_vals[simi_vals>threshold]
                      *dist_weights[simi_vals>threshold])
  }
  result = c()
  result$uncertainty = uncertainty
  result$predicted = predicted
  return(result)
}

