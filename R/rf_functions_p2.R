#' @export
#' @title Smooth AUC Values
#' @description  \code{auc.smooth} takes in a set of time varying AUC values and smooths them out, in addition to providing variance values
#' @param data The dataframe with the target values
#' @param auc The dataframe with the initial AUC values
#' @param target String variable containing the name of the target
#' @param patient_count_col String variable containing the name of the column couting the patient number
#' @return a dataframe with a row for each CPIU and its smoothed AUC value
auc.smooth <- function(data, auc, target, patient_count_col){
  data[, target] <- as.numeric(as.character(data[, target]))
  names(data)[names(data) == patient_count_col] <- "int.n"
  names(data)[names(data) == target] <- "target"
  boot.df.rep <- data

  events.df <- boot.df.rep %>% dplyr::group_by(int.n) %>% dplyr::summarise(pos = sum(target))
  events.df <- as.data.frame(events.df)
  events.df.pos <- events.df %>% filter(pos > 0 )
  auc.df <- auc
  names(auc.df)[1] <- "int.n"
  events.df.pos <- left_join(auc.df, events.df.pos, by = "int.n")
  events.df.pos$auc <- ifelse(events.df.pos$auc>0.95, 0.95, events.df.pos$auc)
  events.df.pos$auc <- ifelse(events.df.pos$auc<0.05, 0.05, events.df.pos$auc)
  events.df.pos <- events.df.pos %>% mutate(var = (auc)*(1-auc)/pos)
  events.df.pos <- events.df.pos %>% mutate(w = 1/var)
  return(events.df.pos)
}

#' @title Weighted Smooth AUC Value
#' @description \code{auc.smooth.return.single} both smooths the AUC values and then provides a single weighted average AUC based on the number
#' of observations used for the AUC value at each time
#' @param data The dataframe with the target values
#' @param auc The dataframe with the initial AUC values
#' @param target String variable containing the name of the target
#' @param patient_count_col String variable containing the name of the column couting the patient number
#' @return a weighted average of the time varying CPIU values using weighting based on the number of patients available in each CPIU
#' @export
auc.smooth.return.single <- function(data, auc, target, patient_count_col){
  data[, target] <- as.numeric(as.character(data[, target]))
  names(data)[names(data) == patient_count_col] <- "int.n"
  names(data)[names(data) == target] <- "target"
  boot.df.rep <- data

  events.df <- boot.df.rep %>% dplyr::group_by(int.n) %>% dplyr::summarise(pos = sum(target))
  events.df <- as.data.frame(events.df)
  events.df.pos <- events.df %>% filter(pos > 0 )
  auc.df <- auc
  names(auc.df)[1] <- "int.n"
  events.df.pos <- left_join(auc.df, events.df.pos, by = "int.n")
  events.df.pos$auc <- ifelse(events.df.pos$auc>0.95, 0.95, events.df.pos$auc)
  events.df.pos$auc <- ifelse(events.df.pos$auc<0.05, 0.05, events.df.pos$auc)
  events.df.pos <- events.df.pos %>% mutate(var = (auc)*(1-auc)/pos)
  events.df.pos <- events.df.pos %>% mutate(w = 1/var)
  weighted_avg_auc <- sum((events.df.pos$auc * events.df.pos$num_individuals) / sum(events.df.pos$num_individuals))
  return(weighted_avg_auc)
}

#' @title Calculate time varying AUC dataframe
#' @description \code{rf.auc} calculates the time varying AUC values and returns them in a dataframe
#' @param sca1.df the data
#' @param target String variable containing the name of the target
#' @param patient_count_col String variable containing the name of the column couting the patient number
#' @param time_col String variable containing the name of the column holding the CPIU number
#' @return time varying AUC dataframe with one column for the CPIU and another for the AUC
#' @export
rf.auc <- function(sca1.df = data, target, patient_count_col, time_col){
  names(sca1.df)[names(sca1.df) == patient_count_col] <- "int.n"
  names(sca1.df)[names(sca1.df) == target] <- "target"
  names(sca1.df)[names(sca1.df) == time_col] <- "q6"
  index <- which(sca1.df$target == 1)# find all intervals where sca occurs
  times <- sca1.df[index, "q6"] # obtain times
  status <- sca1.df[,"target"] # 0/1 indicators
  int <- sca1.df[index, "int.n"] # interval numbers for sca
  n <- length(index)

  int <- int[order(times)] # order the event times
  int <- unique(int)
  n <- length(int)

  auc.df.p <- data.frame(time = int, auc = rep(NA, n), num_individuals = rep(NA, n))

  for(i in 1:n){

    int.index <- which(sca1.df$int.n == int[i]) # only consider individuals at risk at the current interval in consideration
    chf <- sca1.df[int.index, c("pid", "target", "p.hat")]
    auc.df.p[i,2] <- auc(chf$target, chf$p.hat)
    auc.df.p[i,3] <- nrow(chf)
  }

  return(auc.df.p)
}


#' @title Predicting Event Rates for New Data
#' @description \code{risk.adjust.new.be} calculates risk event rates on a new set of test data
#' @param rf the random forest object
#' @param status the target variable corresponding to "training" CPIUs used to build the rf
#' @param rt the risk time values corresponding to "training" CPIUs used to build the rf
#' @param new_data the new test dataset
#' @return vector of risk values for the new data
#' @export
risk.adjust.new.be <- function(rf, status, rt, new_data, k = 2, alpha.tm = 0.05){

  pre.mem <- rfSLAM::predict.rfsrc(rf, new_data, na.action="na.impute", membership = TRUE) #get the predicted values
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

  t.1 <- data.frame(node = pre.mem$membership[,1]) # tree 1 terminal nodes for new data

  t.1.p <- left_join(t.1, p.df, by = "node") # determine predicted event rate for each new CPIU

  cpius <- dim(pre.mem$membership)[1] # number of new CPIUs

  p.hat.mem <- data.frame(matrix(NA, nrow = cpius, ncol = t))
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

    t.1 <- data.frame(node = pre.mem$membership[,i])
    t.1.p <- left_join(t.1, p.df, by = "node")

    p.hat.mem[,i] <- t.1.p$p.hat
  }


  a.p.hat.mem <- apply(p.hat.mem, 1, mean, trim = alpha.tm, na.rm = TRUE)


  return(a.p.hat.mem) # predicted event rates

}

#' @export
#' @title Create Analysis Plots
#' @description \code{analysis_plots} creates plots of both the calibration of the random forest and the rpart summary tree
#' @param rf.df.1 the dataframe
#' @param target the column name of the target variable
#' @param id_Col the column name of the patient id variable
#' @param time_col the column name with the cpiu values
#' @param vars_list the variables used in training the random forest
#' @return shows the calibration and rpart plots for the model
analysis_plots <- function(rf.df.1, target, id_col, risk_col, time_col, vars_list) {
  set.seed(321)

  db2 <- rf.df.1 %>% mutate(ni.sca = (as.numeric(rf.df.1[,target]) - 1)) %>% mutate(i.sca = rf.df.1[,target]) %>% mutate(timescd = as.numeric(rf.df.1[,time_col])) %>% mutate(int.n = as.numeric(rf.df.1[,time_col]))
  db2$pid <- rf.df.1[,id_col]
  db2$p.hat <- rf.df.1[,risk_col]

  g1 <- mutate(db2, bin = ntile(p.hat, 10)) %>%
    # Bin prediction into 10ths
    group_by(bin) %>%
    dplyr::mutate(n = n(), # Get ests and CIs
                  bin_pred = mean(p.hat),
                  bin_prob = mean(as.numeric(ni.sca)),
                  se = sqrt((bin_prob * (1 - bin_prob)) / n),
                  ul = bin_prob + 1.96 * se,
                  ll = bin_prob - 1.96 * se) %>%
    ungroup() %>%
    ggplot(aes(x = bin_pred, y = bin_prob, ymin = ll, ymax = ul)) +
    geom_pointrange(size = 0.5, color = "black") +
    scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 1, by = 0.1)) +
    scale_x_continuous(limits = c(0, 0.6), breaks = seq(0, 1, by = 0.1)) +
    geom_abline() + geom_smooth(method = "loess", se = FALSE, linetype = "dashed",
                                color = "black") +
    xlab("") +
    ylab("Observed Probability") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5,face = "bold")) +
    ggtitle("Calibration")

  g2 <- ggplot(db2, aes(x = p.hat)) +
    geom_histogram(fill = "black", bins = 200) +
    scale_x_continuous(limits = c(0, 0.6), breaks = seq(0, 1, by = 0.1)) +
    xlab("Predicted Probability") +
    ylab("") +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 25000),breaks = seq(0, 25000, by = 5000)) +
    theme(panel.grid.minor = element_blank())

  g <- arrangeGrob(g1, g2, respect = TRUE, heights = c(1, 0.25), ncol = 1)

  grid.newpage()
  grid.draw(g)


  ##Rpart Plot
  analysis_vec <- vars_list
  tree.df <- db2[,intersect(c(analysis_vec, "p.hat"), colnames(db2))]
  global.tree.1 <- rpart(p.hat ~., data = tree.df, control = rpart.control(cp = 0.005))
  tree.global.predict <- predict(global.tree.1, newdata = tree.df)
  preds <- tree.global.predict
  actual <- tree.df$p.hat
  rpart.plot(global.tree.1, box.palette = "Reds", extra = 1, cex = 0.6, type = 5)
}

#' @export
#' @title Create Feature Importance Plot
#' @description \code{feature_importance_plot} Plots the feature importance for the random forest
#' @param mymodel.1.full the random forest model
#' @param var_key a variable key dataframe linking the variable names and their true meaning
#' @param importance_threshold minimum percent of trees variable should be in in order to be included in the plot
#' @return shows the feature importance plot for the model
feature_importance_plot <- function(mymodel.1.full, var_key, importance_threshold = 10) {
  vars.tree <- mymodel.1.full$var.used
  var.tree.sumC <- colSums(vars.tree)
  var.tree.sumR <- rowSums(vars.tree)
  for(m in 1:nrow(vars.tree)){
    for(n in 1:ncol(vars.tree)){
      vars.tree[m,n] <- ifelse(vars.tree[m,n]>1,1,vars.tree[m,n])
    }
  }
  var.tree.sumC <- colSums(vars.tree)
  var.tree.sumR <- rowSums(vars.tree)
  var.used <- data.frame(Variable = names(var.tree.sumC), perc_tree = unname(var.tree.sumC)) %>% filter(perc_tree>importance_threshold)

  var.used <- left_join(var.used,var_key, by = "Variable")

  p1 <- ggplot(data = var.used, aes(x=reorder(Key,perc_tree), y = perc_tree, fill = Type)) +  geom_bar(stat="identity", color = "black") + theme_bw() + coord_flip()

  p1 + scale_fill_npg() + scale_y_continuous(limits = c(0,100), expand = c(0,0), breaks = seq(0, 100, by = 10)) + xlab("") +
    ylab("Percent of Trees Using Variable") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5,face = "bold"), axis.text=element_text(size=6)) +
    ggtitle("Variable Importance")
}


#' @export
#' @title Train an RFSLAM model
#' @description \code{create_model} trains a RFSLAM model using the given user parameter values
#' @param modeling_df dataframe containing data for modeling
#' @param target name of the target variable column
#' @param id_col name of the column with the patient id's
#' @param risk_time_col name of the column with the risk time values
#' @param patient_count_col name of the column with the patient counts
#' @param ntree number of trees for random forest
#' @param nodedepth the node depth for random forest
#' @param nsplit the nsplit parameter for random forest
#' @return the RF-SLAM model, the predicted risk values (uncalibrated) together in a list that you can index into using the `model` and `preds` arguments
create_model <- function(modeling_df, target, id_col, risk_time_col, patient_count_col = "int.n", time_col, ntree = 100, nodedepth = NULL, nsplit = 10) {

  rfslam.samp.1 <- boot.samp(ntree = ntree, id = unique(modeling_df[,id_col]), boot.df.rep = modeling_df %>% select(-c(risk_time_col)), id_col = id_col, time_unit_col = patient_count_col)
  node.size <- round((dim(modeling_df %>% select(-c(id_col, risk_time_col)))[1])*0.1)

  mymodel.1.full <- rfSLAM::rfsrc(as.formula(paste(target, "~", ".")), data = modeling_df
                                  %>% select(-c(id_col, risk_time_col, patient_count_col)), nodesize = node.size, ntree =  ntree, nodedepth = nodedepth, nsplit = nsplit,
                                  na.action = "na.impute", splitrule= "poisson.split1", risk.time = modeling_df[, risk_time_col],
                                  stratifier=modeling_df[,time_col], membership=TRUE,bootstrap = "by.user", samp = rfslam.samp.1, var.used = "by.tree", importance = TRUE)

  p.cpiu.be.ppl <- risk.adjust.be(rf = mymodel.1.full, status = modeling_df[,target], rt = modeling_df[,risk_time_col], alpha.tm = 0)

  return(list("model" = mymodel.1.full, "preds" = p.cpiu.be.ppl))
}

#' @export
#' @title Perform Cross Validation for RFSLAM model
#' @description \code{cv_model_with_auc} performs cross validation and returns the average auc value across the folds for the model
#' @param modeling_df dataframe containing data for modeling
#' @param target name of the target variable column
#' @param id_col name of the column with the patient id's
#' @param risk_time_col name of the column with the risk time values
#' @param patient_count_col name of the column with the patient counts
#' @param n.folds the number of folds for cross validation
#' @param folds_stratifier the variable to use for creating the folds
#' @param drop the variables that need to be dropped before training the model
#' @param ntree number of trees for random forest
#' @param nodedepth the node depth for random forest
#' @param nsplit the nsplit parameter for random forest
#' @return the average weighted auc value across all of the folds
cv_model_with_auc <- function(modeling_df, target, id_col, risk_time_col, patient_count_col = "int.n", time_col, n.folds,folds_stratifier, drop, ntree = 100, nodedepth = NULL, nsplit = 10) {
  all_auc <- c()
  data_for_folds <-  modeling_df %>% select(!!folds_stratifier, !!id_col) %>% distinct(.)
  folds <- caret::createFolds(data_for_folds[,folds_stratifier], k = n.folds, list = FALSE)
  data_for_folds$fold <- folds
  modeling_df <- modeling_df %>% left_join(data_for_folds %>% select(-c(!!folds_stratifier)), by = id_col)

  fold_number <- 1:n.folds
  drop <- c(drop, "fold")
  drop_check <<- drop
  get_fold_auc <- function(j, .modeling_df = modeling_df, .target = target, .id_col = id_col, .drop = drop, .risk_time_col = risk_time_col, .time_col = time_col, .patient_count_col = patient_count_col, .ntree = ntree, .nodedepth = nodedepth, .nsplit = nsplit) {
    train <- .modeling_df %>% filter(fold == j)
    test <- .modeling_df %>% filter(fold != j)
    test <- rbind(train[1, ] , test)
    test <- test[-1,]

    model_pieces <- create_model(train[,!(names(train) %in% .drop[.drop != sym(.target)])], .target, .id_col, .risk_time_col , .patient_count_col, .time_col, .ntree, .nodedepth, .nsplit)
    rf <- model_pieces$model
    test <-  test %>% select(-c(fold))


    # train_check <<- train
    # forest <<- rf
    # test_check <<- test
    # drop_check <<- .drop

    test$p.hat <- risk.adjust.new.be(rf, train[,.target], train[,.risk_time_col], test[,!(names(test) %in% c(.drop, .target))])
    test_check <<- test
    auc_fold <- auc.smooth.return.single(test, rf.auc(test, .target, .patient_count_col, .time_col), .target, .patient_count_col)
    return(auc_fold)
  }


  all_auc <- map(.x = fold_number, .f = get_fold_auc)
  return(all_auc)
}


#this has now been optimized using map
#' @export
#' @title Tune Parameters for RFSLAM
#' @description \code{tune_rf_params} reports model auc values for parameter values searched over a grid as defined by the passed in parameters
#' @param df the dataframe for modeling
#' @param target name of the target variable column
#' @param id_col name of the column with the patient id's
#' @param risk_time_col name of the column with the risk time values
#' @param patient_count_col name of the column with the patient counts
#' @param drop the variables to drop from the dataframe before training the model
#' @param ntree_trys the ntree options to try
#' @param nodedepth_trys the nodedepth options to try
#' @param nsplit_trys the nsplit options to try
#' @param n.folds the number of folds to use for cross validation
#' @param folds_stratifier the variable to use for creating the folds
#' @return a table showing all combinations of the parameter values and their associated cross validated weighted auc value
tune_rf_params <- function(df, target, id_col, risk_time_col, patient_count_col = "int.n", time_col, drop, ntree_trys = c(50, 100, 200, 500), nodedepth_trys = c(NULL, 3, 5), nsplit_trys = c(5, 10, 15), n.folds = 5, folds_stratifier) {

  parameter_combos <- expand.grid(n_trees = ntree_trys, nodedepth = nodedepth_trys, nsplit = nsplit_trys)

  get_combo_auc <- function(n_tree, nodedepth, nsplit, .df = df, .target = target, .id_col = id_col, .risk_time_col = risk_time_col,
                            .patient_count_col = patient_count_col, .time_col = time_col, .n.folds = n.folds, .drop = drop, .folds_stratifier = folds_stratifier) {
    if (nodedepth == "NULL") {
      nodedepth <- NULL
    }
    else {
      nodedepth <- as.numeric(as.character(nodedepth))
    }
    curr_model_auc <- cv_model_with_auc(.df, .target, .id_col, .risk_time_col, .patient_count_col, .time_col, n.folds = .n.folds,
                                        drop = .drop, folds_stratifier = .folds_stratifier, ntree = n_tree, nodedepth = nodedepth, nsplit = nsplit)

    avg_curr_auc <- mean(unlist(curr_model_auc))
    return(avg_curr_auc)
  }


  all_auc <- pmap(list(parameter_combos$n_trees, parameter_combos$nodedepth, parameter_combos$nsplit), .f = get_combo_auc)
  parameter_combos[,"average auc across folds"] <- unlist(all_auc)
  return(parameter_combos)

}

#' @export
#' @title Calculate Risk Times
#' @description \code{calc_risk_times} calculates the risk times for each patient based on the given window lengths
#' @param df the input dataframe
#' @param end_time_col column containing the final time value a patient had data for
#' @param curr_time_col column containing the time value that the row of data was collected
#' @param window_length lenghts of the CPIU intervals in the model
#' @param window_normalizer value to divide by
#' @param col_name name you want the new column created to have
calc_risk_times <- function(df, end_time_col, curr_time_col, window_length, window_normalizer, col_name) { #comment this out

  df[,col_name] <- ((df[,end_time_col] + 1) - df[,curr_time_col]) * window_length / window_normalizer
  bool_col <- paste(col_name, "p", sep = "")
  df[,bool_col] <- ifelse(df %>% pull(col_name) > 1, 1, df %>% pull(col_name))
  return(df)

}


#' @title Save Risk Values
#' @description \code{save_risk_vals} saves the risk values in \code{df} to the desired \code{file_path}
#' @param df the dataframe with the risk values
#' @file_path path to save to
#' @export
save_risk_vals <- function(df, file_path) {
  saveRDS(df, file_path)
}


