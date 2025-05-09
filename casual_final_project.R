rm(list=ls())
library(dplyr)
library(pheatmap)
library(tidyr)
library(ggplot2)
library(stringr)



##################### Data Pre #######################
pet_img_data <- read.csv("UCBERKELEY_AMY_6MM_07Apr2025.csv")
pet_img_data <-pet_img_data %>%
  group_by(PTID) %>%
  filter(SCANDATE == min(SCANDATE)) %>%
  ungroup()
selected_cols <- grep("^CTX.*SUVR$", colnames(pet_img_data), value = TRUE)
pet_img_data <- pet_img_data [, c(colnames(pet_img_data)[1:10], selected_cols)]
colnames(pet_img_data)[-5] <- paste0("pet_", colnames(pet_img_data)[-5])



tau_img_data <- read.csv("UCBERKELEY_TAU_6MM_10Apr2025.csv")
tau_img_data<-tau_img_data  %>%
  group_by(PTID) %>%
  filter(SCANDATE == min(SCANDATE)) %>%
  ungroup()
selected_cols <- grep("^CTX.*SUVR$", colnames(tau_img_data), value = TRUE)
tau_img_data<- tau_img_data[, c(colnames(tau_img_data)[1:10], selected_cols)]
colnames(tau_img_data)[-5] <- paste0("tau_", colnames(tau_img_data)[-5])




demo_data<-read.csv("PTDEMOG_07Apr2025.csv")
demo_data <- demo_data %>%
  dplyr::select(PTID, PTGENDER, PTDOBYY,PTETHCAT,PTEDUCAT)%>%
  dplyr::filter(PTGENDER != -4,
                PTETHCAT != -4 & PTETHCAT != 3) %>%
  dplyr::mutate(PTGENDER = as.factor(PTGENDER),
                PTETHCAT = as.factor(PTETHCAT))


apoe_data<-read.csv("APOERES_10Apr2025.csv")%>%
  dplyr::select(PTID, GENOTYPE)%>%
  dplyr::mutate(GENOTYPE = as.factor(GENOTYPE))%>%
  mutate(apoe_group = case_when(
    str_count(GENOTYPE, "4") == 2 ~ "apoe42",  # two 4's
    str_count(GENOTYPE, "4") == 1 ~ "apoe41",   # one 4
    TRUE ~ "no4"                          # no 4
  ))

manu_data<-read.csv("All_Subjects_PET_Images_18Apr2025.csv")%>%
  dplyr::select(subject_id, pet_mfr)%>%
  dplyr::mutate(scanner = pet_mfr,
                PTID =  subject_id)%>%
  dplyr::select(PTID, scanner)
manu_data<-manu_data %>%
  distinct(PTID, .keep_all = TRUE)

diag_data<-read.csv("DXSUM_27Apr2025.csv")%>%
  dplyr::select(PTID,DIAGNOSIS)%>%
  distinct(PTID, .keep_all = TRUE)

diag_data$DIAGNOSIS<-factor(diag_data$DIAGNOSIS,
       levels = c(1, 2, 3),
       labels = c("CN", "MCI", "Dementia"))

final_data <- demo_data  %>%
  left_join(pet_img_data, by = "PTID") %>%
  mutate(
    pet_SCANDATE = as.Date(pet_SCANDATE),
    age = as.numeric(format(pet_SCANDATE, "%Y")) - PTDOBYY
  ) %>%
  left_join(tau_img_data, by = "PTID") %>%
  left_join(apoe_data, by = "PTID") %>%
  left_join(manu_data, by = "PTID") %>%
  left_join(diag_data, by = "PTID") %>%
  filter(!is.na(pet_SITEID))

id.na <- which(apply(final_data, 1, function(x) any(is.na(x))))
final_data<- final_data[-id.na, ]

write.csv(final_data, "Casual_final_data.csv", row.names = FALSE)

final_data <- read.csv("Casual_final_data.csv")

## scanner type
unique(final_data$scanner)
table(final_data$scanner)

final_data <- final_data %>%
  mutate(scanner = case_when(
    scanner %in% c("GE MEDICAL SYSTEMS", "GEMS") ~ "GEMS",
    scanner %in% c("Philips Medical Systems", "Philips") ~ "Philips",
    scanner %in% c("SIEMENS", "Siemens ECAT", "Siemens/CTI") ~ "Siemens",
    TRUE ~ scanner  # keep original if no match
  ))

## consider whether delete CPS N=17 

final_data$AMY <- rowMeans(final_data[,15:116 ], na.rm = TRUE)
final_data$TAU <- rowMeans(final_data[, 127:228], na.rm = TRUE)

final_data <- final_data %>%
  mutate(
    APOE41 = ifelse(apoe_group == "apoe41", 1, 0),
    APOE42 = ifelse(apoe_group == "apoe42", 1, 0)
  )

final_data <- final_data %>%
  mutate(
    DIAGNOSIS = case_when(
      DIAGNOSIS == "CN" ~ "CN",
      DIAGNOSIS %in% c("MCI", "Dementia") ~ "AD"
    ),
    DIAGNOSIS = factor(DIAGNOSIS, levels = c("CN", "AD"))
  )

write.csv(final_data, "new_Casual_final_data.csv", row.names = FALSE)

##################### FCI model #######################


library(pcalg)
library(dplyr)
library(tibble)
library(Matrix)
library(parallel)

final_data <- read.csv("new_Casual_final_data.csv")

#final_data$old_DIAGNOSIS<-final_data$DIAGNOSIS
final_data$DIAGNOSIS<- ifelse(final_data$DIAGNOSIS=="AD",1,0)
# Select columns and clean
df_clean <- final_data %>%
  dplyr::select(PTGENDER, PTEDUCAT, age,PTETHCAT, APOE41, APOE42, TAU, AMY, DIAGNOSIS) %>%
  mutate(across(everything(), as.numeric)) %>%
  na.omit()

df_clean$PTGENDER<-df_clean$PTGENDER-1
df_clean$PTETHCAT<-df_clean$PTETHCAT-1

var_names <- colnames(df_clean)
num_vars <- length(var_names)
n <- nrow(df_clean)

# Function to apply FCI to one bootstrap sample
run_fci <- function(data) {
  boot_sample <- data[sample(1:nrow(data), replace = TRUE), ]
  suffStat <- list(C = cor(boot_sample), n = nrow(boot_sample))
  fci_fit <- fci(suffStat, indepTest = gaussCItest, alpha = 0.05,
                 labels = var_names, verbose = FALSE)
  return(fci_fit@amat)  # Return adjacency matrix
}

# Run in parallel
set.seed(123)
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("df_clean", "var_names", "fci", "gaussCItest", "run_fci"))
clusterEvalQ(cl, library(pcalg))

fci_results <- parLapply(cl, 1:1000, function(i) run_fci(df_clean))
stopCluster(cl)


# Initialize matrix to count edges
edge_freq <- matrix(0, nrow = num_vars, ncol = num_vars,
                    dimnames = list(var_names, var_names))

# Count presence of each directed edge
for (amat in fci_results) {
  edge_freq <- edge_freq + (amat != 0)
}

# Normalize by number of bootstraps
edge_prob <- edge_freq / length(fci_results)

threshold <- 0.85

cat("Edges appearing in >", threshold * 1000, "% of bootstrap samples:\n")
for (i in 1:num_vars) {
  for (j in 1:num_vars) {
    if (edge_prob[i, j] >= threshold) {
      cat(var_names[i], "→", var_names[j], "with probability", round(edge_prob[i, j], 2), "\n")
    }
  }
}


library(ggplot2)
library(reshape2)

melted <- melt(edge_prob)
ggplot(melted, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Edge Probability Across 100 Bootstrap FCIs",
       x = "To", y = "From") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##################### FCI model evaluation #######################

data<-final_data
#data$APOE4 <- as.integer(data$APOE41 == 1 | data$APOE42 == 1)
data$PTGENDER <- as.factor(data$PTGENDER)
data$DIAGNOSIS <- as.factor(data$DIAGNOSIS)




set.seed(123)
idx <- sample(seq_len(nrow(data)), size = 0.7 * nrow(data))
train <- data[idx, ]
test <- data[-idx, ]

model_bin <- glm(DIAGNOSIS ~ AMY + TAU + PTGENDER + APOE41, 
                 data = train, family = binomial)
summary(model_bin)


# For binary
pred_probs_bin <- predict(model_bin, newdata = test, type = "response")
pred_labels_bin <- ifelse(pred_probs_bin > 0.5, 1, 0)


library(pROC)
roc_obj <- roc(test$DIAGNOSIS, pred_probs_bin)
plot(roc_obj)
auc(roc_obj)


##################### Evaluation intra-site #######################
##Note: last scanner with small sample size, so remove it so far
intra_results <- list()

for (site_i in unique(final_data$scanner)[1:3]) {
  
  cat("Running intra-site for scanner:", site_i, "\n")
  
  data <- final_data %>%
    dplyr::select(PTGENDER, PTEDUCAT, age,PTETHCAT, APOE41, APOE42, TAU, AMY, DIAGNOSIS) %>%
    mutate(across(everything(), as.numeric)) %>%
    na.omit()
  
  data$PTGENDER<-data$PTGENDER-1
  data$PTETHCAT<-data$PTETHCAT-1
  
  data$scanner<- final_data$scanner
  
  site_data <- subset(data, scanner == site_i)
  
  set.seed(123)
  idx <- sample(1:nrow(site_data), 0.7 * nrow(site_data))
  train_data <- site_data[idx, ]
  test_data <- site_data[-idx, ]
  
  # FCI on train data
  suffStat <- list(C = cor(train_data[, c("PTGENDER", "PTEDUCAT", "age","PTETHCAT", "APOE41", "APOE42", "TAU", "AMY", "DIAGNOSIS")], use = "pairwise.complete.obs"),
                   n = nrow(train_data))
  fci_fit <- fci(suffStat, indepTest = gaussCItest, alpha = 0.05, 
                 labels = colnames(suffStat$C), skel.method = "stable")
  
  # Select predictors based on FCI results (edges into DIAGNOSIS_BIN)
  diagnosis_parents <- fci_fit@amat[, "DIAGNOSIS"]
  predictors <- names(which(diagnosis_parents != 0))
  
  # Fit logistic model using predictors from FCI
  if (!length(predictors) == 0){
    train_data$PTGENDER <- as.factor(train_data$PTGENDER)
    train_data$DIAGNOSIS <- as.factor(train_data$DIAGNOSIS)
    
    test_data$PTGENDER <- as.factor(test_data$PTGENDER)
    test_data$DIAGNOSIS <- as.factor(test_data$DIAGNOSIS)
    
    formula_str <- paste("DIAGNOSIS ~", paste(predictors, collapse = " + "))
    model <- glm(as.formula(formula_str), data = train_data, family = binomial)
    
    # Predict on test set
    probs <- predict(model, newdata = test_data, type = "response")
    roc_obj <- roc(test_data$DIAGNOSIS, probs)
    
    intra_results[[site_i]] <- list(
      site = site_i,
      predictors = predictors,
      auc = auc(roc_obj)
    )
  }else{
    print("there is no causal model for this site")
  }
}

##################### Evaluation inter-site #######################

inter_results <- list()
sites <- unique(final_data$scanner)[1:3]

for (site_train in sites) {
  for (site_test in setdiff(sites, site_train)) {
    cat("Train on site:", site_train, "| Test on site:", site_test, "\n")
    
    data <- final_data %>%
      dplyr::select(PTGENDER, PTEDUCAT, age,PTETHCAT, APOE41, APOE42, TAU, AMY, DIAGNOSIS) %>%
      mutate(across(everything(), as.numeric)) %>%
      na.omit()
    
    data$scanner<- final_data$scanner
    
    data$PTGENDER<-data$PTGENDER-1
    data$PTETHCAT<-data$PTETHCAT-1
    
    train_data <- subset(data, scanner == site_train)
    test_data <- subset(data, scanner == site_test)
    
    # FCI on train data
    suffStat <- list(C = cor(train_data[, c("PTGENDER", "PTEDUCAT", "age","PTETHCAT", "APOE41", "APOE42", "TAU", "AMY", "DIAGNOSIS")], use = "pairwise.complete.obs"),
                     n = nrow(train_data))
    fci_fit <- fci(suffStat, indepTest = gaussCItest, alpha = 0.05,
                   labels = colnames(suffStat$C), skel.method = "stable")
    
    diagnosis_parents <- fci_fit@amat[, "DIAGNOSIS"]
    predictors <- names(which(diagnosis_parents != 0))
    
    if (!length(predictors) == 0){
    formula_str <- paste("DIAGNOSIS ~", paste(predictors, collapse = " + "))
    model <- glm(as.formula(formula_str), data = train_data, family = binomial)
    
    probs <- predict(model, newdata = test_data, type = "response")
    roc_obj <- roc(test_data$DIAGNOSIS, probs)
    
    inter_results[[paste(site_train, "->", site_test)]] <- list(
      train_site = site_train,
      test_site = site_test,
      predictors = predictors,
      auc = auc(roc_obj)
    )
    }else{
      print("there is no causal model")
    }
  }}



##################### LiNGAM model #######################

library(glmnet)
library(stabs)


# Step 0: Convert DIAGNOSIS to numeric
X_full<- final_data[,c(2,4,5,15:116,127:228,232,235,236)]


X_full$DIAGNOSIS <- as.factor(X_full$DIAGNOSIS)
X_full$PTGENDER <-X_full$PTGENDER-1 
X_full$PTEDUCAT <-X_full$PTEDUCAT-1 

# Separate predictors and outcome
y <- X_full$DIAGNOSIS
X <- X_full[, setdiff(colnames(X_full), "DIAGNOSIS")]


set.seed(123)

# Selected variable indices and names

lambda.min <- cv.glmnet(x = as.matrix(X), y = y, family = "binomial")$lambda.min
stab.maxCoef <- stabsel(x = X,
                         y = y,
                         fitfun = glmnet.lasso_maxCoef, 
                         # specify additional parameters to fitfun
                         args.fitfun = list(lambda = lambda.min),
                         cutoff = 0.75, PFER = 1)

selected_vars <- colnames(X)[stab.maxCoef$selected]
print(selected_vars)

X_sub<- X[, selected_vars]
#X_sub$PTGENDER<-final_data$PTGENDER-1
#X_sub$PTETHCAT<-final_data$PTETHCAT-1
#X_sub$PTEDUCAT<-final_data$PTEDUCAT
#X_sub$APOE41<-final_data$APOE41
#X_sub$APOE42<-final_data$APOE42
X_sub$DIAGNOSIS<-X_full$DIAGNOSIS

X_sub <- X_sub %>%mutate(across(everything(), as.numeric)) %>%
  na.omit()

# Step 4: Apply LiNGAM to connected variables
#library(Rlingam)
lingam.result <- lingam(X_sub)

# Step 5: View the estimated causal adjacency matrix
print("LiNGAM Estimated Causal Matrix:")
print(lingam.result$B)

# Optional: visualize
library(igraph)

B<-lingam.result$B
colnames(B)<- rownames(B)<-colnames(X_sub)

## plot function have bug
plot_lingam_graph <- function(B_matrix, threshold = 1e-6, title = "LiNGAM Causal Graph") {
  if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
  library(igraph)
  
  # Check for row/col names
  if (is.null(colnames(B_matrix)) || is.null(rownames(B_matrix))) {
    stop("Please provide row and column names for B_matrix.")
  }
  
  # Extract nonzero edges (above threshold)
  edges <- which(abs(B_matrix) > threshold, arr.ind = TRUE)
  
  if (nrow(edges) == 0) {
    warning("No edges found in B_matrix above threshold.")
    return(NULL)
  }
  
  # Create edge list
  edge_list <- data.frame(
    from   = colnames(B_matrix)[edges[, 2]],
    to     = rownames(B_matrix)[edges[, 1]],
    weight = B_matrix[edges]
  )
  
  all_nodes <- colnames(B_matrix)
  
  
  g <- graph_from_data_frame(
    d = edge_list,
    directed = TRUE,
    vertices = data.frame(name = all_nodes)  # Include all nodes
  )
  
  # Assign edge color by sign and add labels
  E(g)$color <- ifelse(edge_list$weight > 0, "darkgreen", "red")
  E(g)$label <- round(edge_list$weight, 3)  # Edge weights as labels
  
  # Layout using all nodes (even disconnected ones)
  layout <- layout_with_fr(
    graph = g,
    weights = abs(E(g)$weight)  # Use absolute weights for layout
  )
  
  # Plot
  plot(g,
       layout = layout,
       vertex.size = 35,
       vertex.color = "skyblue",
       vertex.label = V(g)$name,  # Use node names from graph
       vertex.label.cex = 1.2,
       vertex.label.color = "black",
       edge.arrow.size = 0.8,
       edge.width = 2,
       edge.label.cex = 1.1,      # Larger edge labels
       edge.label.color = "black",
       main = title)
}

plot_lingam_graph(B)

##################### LiNGAM model evaluation #######################

data<-final_data
#data$APOE4 <- as.integer(data$APOE41 == 1 | data$APOE42 == 1)
data$PTGENDER <- as.factor(data$PTGENDER)
data$DIAGNOSIS <- as.factor(data$DIAGNOSIS)




set.seed(123)
idx <- sample(seq_len(nrow(data)), size = 0.7 * nrow(data))
train <- data[idx, ]
test <- data[-idx, ]

model_bin <- glm(DIAGNOSIS ~ PTGENDER + pet_CTX_RH_ENTORHINAL_SUVR + PTGENDER + tau_CTX_LH_ROSTRALANTERIORCINGULATE_SUVR, 
                 data = train, family = binomial)
summary(model_bin)


# For binary
pred_probs_bin <- predict(model_bin, newdata = test, type = "response")
pred_labels_bin <- ifelse(pred_probs_bin > 0.5, 1, 0)


library(pROC)
roc_obj <- roc(test$DIAGNOSIS, pred_probs_bin)
plot(roc_obj)
auc(roc_obj)


##################### Evaluation intra-site #######################
##Note: last scanner with small sample size, so remove it so far
intra_results <- list()

for (site_i in unique(final_data$scanner)[1:3]) {
  
  cat("Running intra-site for scanner:", site_i, "\n")
  
  data <- final_data
  site_data <- subset(data, scanner == site_i)
  
  set.seed(123)
  idx <- sample(1:nrow(site_data), 0.7 * nrow(site_data))
  train_data <- site_data[idx, ]
  test_data <- site_data[-idx, ]
  
  X_full<- train_data[,c(2,4,5,15:116,127:228,232,235,236)]
  
  
  X_full$DIAGNOSIS <- as.factor(X_full$DIAGNOSIS)
  X_full$PTGENDER <-X_full$PTGENDER-1 
  X_full$PTEDUCAT <-X_full$PTEDUCAT-1 
  
  # Separate predictors and outcome
  y <- X_full$DIAGNOSIS
  X <- X_full[, setdiff(colnames(X_full), "DIAGNOSIS")]
  
  
  set.seed(123)
  
  # Selected variable indices and names
  
  lambda.min <- cv.glmnet(x = as.matrix(X), y = y, family = "binomial")$lambda.min
  stab.maxCoef <- stabsel(x = X,
                          y = y,
                          fitfun = glmnet.lasso_maxCoef, 
                          # specify additional parameters to fitfun
                          args.fitfun = list(lambda = lambda.min),
                          cutoff = 0.75, PFER = 1)
  
  selected_vars <- colnames(X)[stab.maxCoef$selected]
  print(selected_vars)
  
  X_sub<- X[, selected_vars]
  predictors<-colnames(X_sub)
  X_sub$DIAGNOSIS<-X_full$DIAGNOSIS
  
  X_sub <- X_sub %>%mutate(across(everything(), as.numeric)) %>%
    na.omit()
  
  lingam.result <- lingam(X_sub)
  
  print("LiNGAM Estimated Causal Matrix:")
  print(lingam.result$B)
  
  nonzero_cols <- which(lingam.result$B[nrow(lingam.result$B), ] != 0)
  predictors<-predictors[nonzero_cols]
  
  # Fit logistic model using predictors from FCI
  if (!length(predictors) == 0){
    train_data$PTGENDER <- as.factor(train_data$PTGENDER)
    train_data$DIAGNOSIS <- as.factor(train_data$DIAGNOSIS)
    
    test_data$PTGENDER <- as.factor(test_data$PTGENDER)
    test_data$DIAGNOSIS <- as.factor(test_data$DIAGNOSIS)
    
    formula_str <- paste("DIAGNOSIS ~", paste(predictors, collapse = " + "))
    model <- glm(as.formula(formula_str), data = train_data, family = binomial)
    
    # Predict on test set
    probs <- predict(model, newdata = test_data, type = "response")
    roc_obj <- roc(test_data$DIAGNOSIS, probs)
    
    intra_results[[site_i]] <- list(
      site = site_i,
      predictors = predictors,
      auc = auc(roc_obj)
    )
  }else{
    print("there is no causal model for this site")
  }
}

##################### Evaluation inter-site #######################

inter_results <- list()
sites <- unique(final_data$scanner)[1:3]

for (site_train in sites) {
  for (site_test in setdiff(sites, site_train)) {
    cat("Train on site:", site_train, "| Test on site:", site_test, "\n")
    
    ata <- final_data
    site_data <- subset(data, scanner == site_i)
    
    set.seed(123)
    idx <- sample(1:nrow(site_data), 0.7 * nrow(site_data))
    train_data <- site_data[idx, ]
    test_data <- site_data[-idx, ]
    
    X_full<- train_data[,c(2,4,5,15:116,127:228,232,235,236)]
    
    
    X_full$DIAGNOSIS <- as.factor(X_full$DIAGNOSIS)
    X_full$PTGENDER <-X_full$PTGENDER-1 
    X_full$PTEDUCAT <-X_full$PTEDUCAT-1 
    
    # Separate predictors and outcome
    y <- X_full$DIAGNOSIS
    X <- X_full[, setdiff(colnames(X_full), "DIAGNOSIS")]
    
    
    set.seed(123)
    
    # Selected variable indices and names
    
    lambda.min <- cv.glmnet(x = as.matrix(X), y = y, family = "binomial")$lambda.min
    stab.maxCoef <- stabsel(x = X,
                            y = y,
                            fitfun = glmnet.lasso_maxCoef, 
                            # specify additional parameters to fitfun
                            args.fitfun = list(lambda = lambda.min),
                            cutoff = 0.75, PFER = 1)
    
    selected_vars <- colnames(X)[stab.maxCoef$selected]
    print(selected_vars)
    
    X_sub<- X[, selected_vars]
    predictors<-colnames(X_sub)
    X_sub$DIAGNOSIS<-X_full$DIAGNOSIS
    
    X_sub <- X_sub %>%mutate(across(everything(), as.numeric)) %>%
      na.omit()
    
    lingam.result <- lingam(X_sub)
    
    print("LiNGAM Estimated Causal Matrix:")
    print(lingam.result$B)
    
    nonzero_cols <- which(lingam.result$B[nrow(lingam.result$B), ] != 0)
    predictors<-predictors[nonzero_cols]
    
    if (!length(predictors) == 0){
      formula_str <- paste("DIAGNOSIS ~", paste(predictors, collapse = " + "))
      model <- glm(as.formula(formula_str), data = train_data, family = binomial)
      
      probs <- predict(model, newdata = test_data, type = "response")
      roc_obj <- roc(test_data$DIAGNOSIS, probs)
      
      inter_results[[paste(site_train, "->", site_test)]] <- list(
        train_site = site_train,
        test_site = site_test,
        predictors = predictors,
        auc = auc(roc_obj)
      )
    }else{
      print("there is no causal model")
    }
  }}




