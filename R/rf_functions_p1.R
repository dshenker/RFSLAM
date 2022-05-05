#' @title Bootstrap Sampling for RFSLAM
#' @description \code{boot.samp} implements bootstrapping for RFSLAM
#' @param ntree the number of trees to sample for
#' @param id the list of patient ids
#' @param boot.df.rep the dataframe to bootstrap from
#' @param id_col the variable name holding the patient ids
#' @param time_unit_col the variable name holding the CPIU value
boot.samp <- function(ntree = 100, id, boot.df.rep, id_col, time_unit_col){
  set.seed(321)
  # for rfslam, want to get array with number of times each CPIU is in bag
  cpius <- dim(boot.df.rep)[1]
  samp <- array(NA,dim=c(cpius, ntree))
  boot.df.rep.cpiu <- boot.df.rep[,c(id_col, time_unit_col)]
  ppl <- length(id)

  for(i in 1:ntree){
    pids <- sample(id, size = ppl, replace = TRUE)
    pid.df <- as.data.frame(table(pids))
    names(pid.df)[1] <- id_col
    boot.df <- left_join(boot.df.rep.cpiu, pid.df, by = id_col)
    boot.df$Freq <- ifelse(is.na(boot.df$Freq), 0, boot.df$Freq)
    samp[,i] <- boot.df$Freq
  }

  return(samp)
}

#' @export
#' @title Carlibrate RFSLAM model
#' @description \code{calibrate.model} calibrates an RFSLAM model using the true event rates
#' @param p.hat the event risk predictions
#' @param rf.df.1 the dataframe used for modeling
#' @param target_varname name of the variable with the target
#' @param time_varname name of the variable with the CPIU count
calibrate.model <- function(p.hat, rf.df.1, target_varname, time_varname) {
  rf.df.1 <- tibble::rowid_to_column(rf.df.1, "ID")
  rf.df.1$p.hat <- p.hat
  rf.df.1$ni.sca <- as.numeric(rf.df.1[,target_varname]) - 1
  rf.df.1$q6 <- rf.df.1[,time_varname]

  shuffled <- rf.df.1[sample(nrow(rf.df.1)),]

  # create 5 equally sized folds
  folds <- cut(seq(1,nrow(shuffled)), breaks = 5, labels = FALSE)
  shuffled$fold <- folds
  shuffled$risk <- NA

  for(i in 1:5){
    testing <- shuffled[shuffled$fold == i,]
    training <- shuffled[shuffled$fold != i,]
    lr_model = glm(as.factor(ni.sca) ~ ns(p.hat,2)*ns(q6,2),data = training, family = binomial)
    p.hat.df <- data.frame(p.hat = testing$p.hat, q6 = testing$q6)
    lr_probs = predict(lr_model,  newdata = p.hat.df, type = "response")

    shuffled[shuffled$fold == i,"risk"] <- lr_probs
  }
  sorted_df <- shuffled[order(shuffled$ID),]
  return(sorted_df$risk)

}


risk.adjust.be <- function(rf, status, rt, k = 2, alpha.tm = 0.05){

  rfc.simple.dat <- rf
  in.mem <- rfc.simple.dat$membership # which terminal node CPIU is in
  in.bag <- rfc.simple.dat$inbag # inbag means CPIU was used to construct the tree
  oob.df <- ifelse(in.bag == 0, 1, 0)
  in.df <- ifelse(oob.df == 0, 1, 0)
  intervals <- dim(in.mem)[1] # number of CPIUs
  t <- dim(in.mem)[2] # number of trees


  # if out of bag, don't use observation in forming estimate
  for(i in 1:intervals){
    for(j in 1:t){
      if(oob.df[i,j] == 1){
        in.mem[i,j] <- NA # NA means not in bag (aka out of bag)
      }
    }
  }

  # status must be numeric
  if(class(status) == "factor"){
    status <- as.numeric(status)-1
  }
  status <- (status - 1)*-1
  alpha = 1/(k^2)
  lambda.hat = sum(status)/sum(rt)
  beta = alpha/lambda.hat

  i <- 1

  mem <- data.frame(event = status, node = in.mem[,i], rt = rt, oob = oob.df[,i], boot.num = in.bag[,i], in.bag = in.df[,i])

  oob.mem <- data.frame(node = rfc.simple.dat$membership[,i], oob = oob.df[,i]) # terminal nodes
  oob.mem$node <- ifelse(oob.mem$oob==1, oob.mem$node, NA)

  term.e <- mem %>% group_by(node) %>% dplyr::summarise(w.events = sum(boot.num*event))
  term.t <- mem %>% group_by(node) %>% dplyr::summarise(total = sum(boot.num*rt))

  # calculate the number of events per node

  term.e <- as.data.frame(term.e)
  term.e <- term.e[!is.na(term.e$node),] # remove NA node (i.e. oob data)

  # calculate total risk time in node
  term.t <- as.data.frame(term.t)
  term.t <- term.t[!is.na(term.t$node),] # remove NA node (i.e. oob data)


  p.df <- left_join(term.e, term.t, by = "node")
  p.df <- p.df %>% mutate(p.hat = (alpha+w.events)/(beta+total))

  # for each CPIU, determine predicted value from terminal node estimate
  t.1 <- data.frame(event = status, node = rfc.simple.dat$membership[,1], oob = oob.df[,1]) # terminal nodes for tree 1
  t.1.p <- left_join(t.1, p.df, by = "node")


  # create p.hat matrix for membership prediction
  p.hat.mem <- data.frame(matrix(NA, nrow = intervals, ncol = t))
  p.hat.mem[,1] <- t.1.p$p.hat

  ntree <- t
  for(i in 2:ntree){

    mem <- data.frame(event = status, node = in.mem[,i], rt = rt, oob = oob.df[,i], boot.num = in.bag[,i], in.bag = in.df[,i])
    oob.mem <- data.frame(node = rfc.simple.dat$membership[,i], oob = oob.df[,i]) # terminal nodes
    oob.mem$node <- ifelse(oob.mem$oob==1, oob.mem$node, NA)

    term.e <- mem %>% group_by(node) %>% dplyr::summarise(w.events = sum(boot.num*event))
    term.t <- mem %>% group_by(node) %>% dplyr::summarise(total = sum(boot.num*rt))

    term.e <- as.data.frame(term.e)
    term.e <- term.e[!is.na(term.e$node),]

    term.t <- as.data.frame(term.t)
    term.t <- term.t[!is.na(term.t$node),]

    p.df <- left_join(term.e, term.t, by = "node")
    p.df <- p.df %>% mutate(p.hat = (alpha+w.events)/(beta+total))

    t.1 <- data.frame(event = status, node = rfc.simple.dat$membership[,i], oob = oob.df[,i]) # terminal nodes for tree i
    t.1.p <- left_join(t.1, p.df, by = "node")
    t.1.p$p.hat <- ifelse(t.1.p$oob==1, t.1.p$p.hat, NA)

    p.hat.mem[,i] <- t.1.p$p.hat
  }


  a.p.hat.mem <- apply(p.hat.mem, 1, mean, trim = alpha.tm, na.rm = TRUE)

  return(a.p.hat.mem) # predicted event rates
}

#' @title Check the Calibration of an RFSLAM model
#' @description \code{check.cal} checks how well calibrated an RFSLAM model is using the predicted event rates and the true events
#' @param predicted.rate the event risk predictions
#' @param actual.outcomes the true events
#' @param rt the risk time values for each patient
check.cal <- function(predicted.rate, actual.outcomes, rt){
  cal.df <- data.frame(p = predicted.rate, o = actual.outcomes, rt = rt)
  names(cal.df) <- c("p", "o", "rt")
  cal.df$d.p <- as.factor(ntile(cal.df$p, 10))
  cal.df.table <- cal.df %>% group_by(d.p) %>% summarise(d.mean = mean(p), events = sum(o), n.obs = length(o), total.rt = sum(rt), f.events = sum(o)/sum(rt))
  cal.df.table <- cal.df.table %>% mutate(expected = d.mean*total.rt)
  return(as.data.frame(cal.df.table))
}
