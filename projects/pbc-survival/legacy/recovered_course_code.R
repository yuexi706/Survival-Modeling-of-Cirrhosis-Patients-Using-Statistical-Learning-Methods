# Recovered from the supplied rendered course HTML.
# Only local data paths were replaced. Original methodological limitations remain.
# Archival reference; NOT the maintained analysis pipeline.
# Team authorship: see ../../../CONTRIBUTIONS.md.

library(tidyverse)

library(survival)
library(survminer)

library(glmnet)

library(randomForestSRC)

library(gbm)

library(cmprsk)
library(ggcorrplot)
library(naniar)
library(forestmodel)
library(viridis)

library(dplyr)


set.seed(2025)

cirr <- read.csv(Sys.getenv("ORIGINAL_PBC_CSV", "data/cirrhosis.csv"))

glimpse(cirr)

cirr_clean <- cirr %>%
  mutate(
    Status  = factor(Status, levels = c("C", "CL", "D")),
    Drug    = factor(Drug),
    Sex     = factor(Sex),
    Ascites = factor(Ascites),
    Hepatomegaly = factor(Hepatomegaly),
    Spiders = factor(Spiders),
    Edema   = factor(Edema),
    Stage   = factor(Stage)
  )

num_vars <- cirr_clean %>% select_if(is.numeric) %>% names()

for (v in num_vars) {
  med_val <- median(cirr_clean[[v]], na.rm = TRUE)
  cirr_clean[[v]][is.na(cirr_clean[[v]])] <- med_val
}

cat_vars <- cirr_clean %>% select_if(where(is.factor)) %>% names()

for (v in cat_vars) {
  tab <- table(cirr_clean[[v]])
  if (length(tab) > 0) {
    mode_val <- names(tab)[which.max(tab)]
    cirr_clean[[v]][is.na(cirr_clean[[v]])] <- mode_val
  }
}

cont_vars <- c("Age","Bilirubin","Cholesterol","Albumin","Copper",
               "Alk_Phos","SGOT","Tryglicerides","Platelets","Prothrombin")

cirr_clean %>%
  .[, cont_vars] %>%
  pivot_longer(cols = everything(), names_to="variable", values_to="value") %>%
  ggplot(aes(x=value))+
  geom_histogram(bins=30,alpha=.7)+
  facet_wrap(~variable,scales="free")+ theme_minimal()

cirr_clean %>%
  pivot_longer(cols = all_of(cont_vars), names_to="variable", values_to="value") %>%
  ggplot(aes(x=Stage,y=value,fill=Stage))+
  geom_boxplot(alpha=.7)+facet_wrap(~variable,scales="free")+
  scale_fill_viridis(discrete=TRUE)+theme_minimal()

corr_mat <- cor(cirr_clean[, cont_vars])
ggcorrplot(corr_mat, hc.order=TRUE, type="lower", lab=TRUE)

vis_miss(cirr)

cirr_surv <- cirr_clean %>%
  mutate(
    event_death = if_else(Status=="D",1,0),
    cr_status = case_when(
      Status=="D" ~ 1L,
      Status=="CL" ~ 2L,
      TRUE ~ 0L
    )
  )

# Stage
png(filename = "KM_Survival_Plot_Stage.png", width = 8, height = 6, units = "in", res = 300)
fit_stage <- survfit(Surv(N_Days, event_death) ~ Stage, data = cirr_surv)

ggsurvplot(
  fit_stage,
  data = cirr_surv,
  pval = TRUE,
  risk.table = TRUE,
  legend.title = "Stage",
  ggtheme = theme_minimal()
)

dev.off()

knitr::include_graphics("KM_Survival_Plot_Stage.png")

# Plot_Sex
png(filename = "KM_Survival_Plot_Sex.png", width = 8, height = 6, units = "in", res = 300)
fit_sex <- survfit(Surv(N_Days, event_death) ~ Sex, data = cirr_surv)

ggsurvplot(
  fit_sex,
  data = cirr_surv,
  pval = TRUE,
  risk.table = TRUE,
  legend.title = "Sex",
  ggtheme = theme_minimal()
)

dev.off()

knitr::include_graphics("KM_Survival_Plot_Sex.png")

# Plot_Drug
png(filename = "KM_Survival_Plot_Drug.png", width = 8, height = 6, units = "in", res = 300)
fit_drug <- survfit(Surv(N_Days, event_death) ~ Drug, data = cirr_surv)

ggsurvplot(
  fit_drug,
  data = cirr_surv,
  pval = TRUE,
  risk.table = TRUE,
  legend.title = "Drug",
  ggtheme = theme_minimal()
)

dev.off()

knitr::include_graphics("KM_Survival_Plot_Drug.png")

n <- nrow(cirr_surv)
train_idx <- sample(seq_len(n),floor(0.7*n))
train_data <- cirr_surv[train_idx,]
test_data  <- cirr_surv[-train_idx,]

surv_train <- with(train_data,Surv(N_Days,event_death))

cox_formula <- as.formula(
  Surv(N_Days,event_death) ~ Drug + Age + Sex + Ascites + Hepatomegaly +
    Spiders + Edema + Bilirubin + Cholesterol + Albumin + Copper + 
    Alk_Phos + SGOT + Tryglicerides + Platelets + Prothrombin + Stage
)

cox_fit <- coxph(cox_formula, data=train_data)
cox_lp_test <- predict(cox_fit,newdata=test_data,type="lp")

cox_cindex <- survConcordance(
  Surv(test_data$N_Days,test_data$event_death) ~ cox_lp_test
)$concordance

step_cox <- step(cox_fit, direction="both", trace=FALSE)
forest_model(step_cox)

cols_to_exclude <- c("N_Days", "Status", "event_death", "cr_status")
x_vars <- cirr_surv[, !(names(cirr_surv) %in% cols_to_exclude)]

x_train <- model.matrix(~.,data=x_vars[train_idx,])[,-1]
x_test  <- model.matrix(~.,data=x_vars[-train_idx,])[,-1]
x_all   <- model.matrix(~.,data=x_vars)[,-1]

cv_lasso <- cv.glmnet(x_train, surv_train, family="cox", alpha=1)
lasso_fit <- glmnet(x_train, surv_train, family="cox", lambda=cv_lasso$lambda.min)

lasso_lp_test <- predict(lasso_fit,newx=x_test,s=cv_lasso$lambda.min)
lasso_cindex <- survConcordance(
  Surv(test_data$N_Days,test_data$event_death) ~ lasso_lp_test
)$concordance

print("--- LASSO-Cox Model: Selected Coefficients at lambda.min ---")

lasso_coefs <- coef(lasso_fit, s = cv_lasso$lambda.min)
print(lasso_coefs[lasso_coefs[, 1] != 0, , drop = FALSE])

# Model C-index
print(paste("LASSO-Cox Test C-index:", round(lasso_cindex, 4)))

plot(cv_lasso)

rsf_fit <- rfsrc(
  Surv(N_Days,event_death) ~ Drug + Age + Sex + Ascites + Hepatomegaly +
    Spiders + Edema + Bilirubin + Cholesterol + Albumin + Copper +
    Alk_Phos + SGOT + Tryglicerides + Platelets + Prothrombin + Stage,
  data=train_data, ntree=1000, importance=TRUE
)

rsf_pred_test <- predict(rsf_fit,newdata=test_data)
rsf_risk_test <- rowMeans(rsf_pred_test$chf)

rsf_cindex <- survConcordance(
  Surv(test_data$N_Days,test_data$event_death) ~ rsf_risk_test
)$concordance

print("--- Random Survival Forest Model Summary ---")

summary(rsf_fit)

# Variable Importance
rsf_vimp <- vimp(rsf_fit)
print(rsf_vimp)

plot(rsf_vimp)

# Model C-index
print(paste("RSF Test C-index:", round(rsf_cindex, 4)))

# Variable Importance (from rfsrc output)
vimp_df <- data.frame(
  Variable = names(rsf_fit$importance),
  VIMP = rsf_fit$importance
)

# Sort by importance
vimp_df <- dplyr::arrange(vimp_df, desc(VIMP))

# Take top 15
vimp_top15 <- vimp_df[1:15, ]

# Plot
ggplot(vimp_top15, aes(x = reorder(Variable, VIMP), y = VIMP)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 15 Variable Importance (Random Survival Forest)",
    x = "Variable",
    y = "VIMP"
  ) +
  theme_minimal()

gbm_fit <- gbm(
  formula = cox_formula,
  data=train_data,
  distribution="coxph",
  n.trees=3000,
  interaction.depth=3,
  shrinkage=0.005,
  cv.folds=5
)

best_trees <- gbm.perf(gbm_fit,method="cv")

gbm_lp_test <- predict(gbm_fit,newdata=test_data,n.trees=best_trees,type="link")

gbm_cindex <- survConcordance(
  Surv(test_data$N_Days,test_data$event_death) ~ gbm_lp_test
)$concordance

# Model Summary
print("--- Gradient-Boosted Model (GBM) Summary ---")

summary(gbm_fit, n.trees=best_trees, cBars = 10)

# Model C-index
print(paste("GBM Test C-index:", round(gbm_cindex, 4)))

ftime <- cirr_surv$N_Days
fstatus <- cirr_surv$cr_status

cov_for_cr <- model.matrix(
  ~ Age + Sex + Bilirubin + Albumin + Prothrombin + Stage,
  data=cirr_surv
)[,-1]

fg_death <- crr(ftime,fstatus,cov1=cov_for_cr,failcode=1,cencode=0)
summary(fg_death)

fg_fit <- cuminc(ftime,fstatus,group=cirr_surv$Stage)
plot(fg_fit)

# 1. Cox Risk Score
cox_lp_all <- predict(cox_fit, newdata = cirr_surv, type = "lp")
cirr_surv$cox_risk <- ntile(cox_lp_all, 3)  # split into tertiles

cox_table <- table(cirr_surv$cox_risk)
print("Cox Risk Stratification Counts:")

print(cox_table)

png(
  filename = "Cox_Risk_Stratification_KM_Curve.png", width = 8, height = 6, units = "in", res = 300)

# Kaplan-Meier plot by Cox risk
cox_surv <- Surv(cirr_surv$N_Days, cirr_surv$cr_status==1)
fit_cox_risk <- survfit(cox_surv ~ cox_risk, data = cirr_surv)
ggsurvplot(fit_cox_risk, data = cirr_surv,
           risk.table = TRUE, pval = TRUE,
           palette = c("steelblue", "orange", "black"),
           title = "Cox Model Risk Stratification")

dev.off()

knitr::include_graphics("Cox_Risk_Stratification_KM_Curve.png")

# 2. LASSO Risk Score

lasso_lp_all <- predict(lasso_fit, newx = x_all, s = cv_lasso$lambda.min)
cirr_surv$lasso_risk <- ntile(lasso_lp_all, 3)

lasso_table <- table(cirr_surv$lasso_risk)
print("LASSO Risk Stratification Counts:")

print(lasso_table)

png(
  filename = "LASSO_Risk_Stratification_KM_Curve.png", width = 8, height = 6, units = "in", res = 300 )
# Kaplan-Meier plot by LASSO risk
fit_lasso_risk <- survfit(cox_surv ~ lasso_risk, data = cirr_surv)
ggsurvplot(fit_lasso_risk, data = cirr_surv,
           risk.table = TRUE, pval = TRUE,
           palette = c("steelblue", "orange", "black"),
           title = "LASSO Model Risk Stratification")

dev.off()

knitr::include_graphics("LASSO_Risk_Stratification_KM_Curve.png")

# 3. RSF Risk Score
rsf_pred_all <- predict(rsf_fit, newdata = cirr_surv)
rsf_lp_all <- rowMeans(rsf_pred_all$chf)  # average cumulative hazard
cirr_surv$rsf_risk <- ntile(rsf_lp_all, 3)

rsf_table <- table(cirr_surv$rsf_risk)
print("RSF Risk Stratification Counts:")

print(rsf_table)

png(
  filename = "RSF_Risk_Stratification_KM_Curve.png", width = 8, height = 6, units = "in", res = 300 )
# Kaplan-Meier plot by RSF risk
fit_rsf_risk <- survfit(cox_surv ~ rsf_risk, data = cirr_surv)
ggsurvplot(fit_rsf_risk, data = cirr_surv,
           risk.table = TRUE, pval = TRUE,
           palette = c("steelblue", "orange", "black"),
           title = "RSF Model Risk Stratification")

dev.off()

knitr::include_graphics("RSF_Risk_Stratification_KM_Curve.png")

cindex_df <- data.frame(
  Model   = c("Cox", "LASSO", "RSF", "GBM"),
  Cindex  = c(cox_cindex, lasso_cindex, rsf_cindex, gbm_cindex)
)

ggplot(cindex_df, aes(x = Model, y = Cindex, fill = Model)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = round(Cindex, 3)), vjust = -0.5) +
  theme_minimal() +
  labs(title = "Model C-index Comparison", y = "C-index") +
  scale_fill_viridis(discrete = TRUE)

risk_df <- data.frame(
  Cox   = cox_lp_all,
  LASSO = as.numeric(lasso_lp_all),
  RSF   = rsf_lp_all
)

risk_long <- pivot_longer(risk_df, cols = everything(), names_to="Model", values_to="Risk")
ggplot(risk_long, aes(x=Model, y=Risk, fill=Model)) +
  geom_violin(alpha=0.7) +
  geom_boxplot(width=0.1, fill="white") +
  theme_minimal(base_size=14) +
  labs(title="Distribution of Risk Scores by Model")

