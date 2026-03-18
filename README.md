# Survival-Modeling-of-Cirrhosis-Patients-Using-Statistical-Learning-Methods

## Overview
### 
This project develops and evaluates machine learning–based survival models to predict mortality risk in patients with liver cirrhosis. The goal is to support clinical decision-making, including patient monitoring, risk stratification, and transplant timing.
Using real-world clinical data from a Mayo Clinic study, we compare traditional statistical approaches with modern machine learning techniques to improve predictive performance and interpretability.

## Dataset
###
- Source: Mayo Clinic cirrhosis study (Primary Biliary Cirrhosis)
- Dataset link (Kaggle): (https://www.kaggle.com/datasets/joebeachcapital/cirrhosis-patient-survival-prediction)
- Sample size: 418 patients
- Features: 17 clinical variables (demographics, symptoms, lab biomarkers)
- Outcome:Death, Censored (alive), Liver transplant (treated as competing risk)

## Methods
###
We implemented and compared multiple survival analysis approaches:

- Cox Proportional Hazards Model
- LASSO-Cox Regression (feature selection + regularization)
- Random Survival Forests (RSF)
- Gradient Boosted Survival Models (GBM)
- Fine-Gray Model (competing risks)

## Key Results
###
- Best Model: LASSO-Cox
  - C-index: 0.857
- All models achieved strong predictive performance (C-index: 0.847–0.857)
- Machine learning models captured nonlinear relationships and improved robustness

## Risk Stratification
###
- We developed a risk scoring system to classify patients into:
  - Low-risk, Medium-risk, High-risk
- Result:
  - Clear survival separation between groups (statistically significant)
  - Enables identification of patients needing:
    - Intensive monitoring
    - Early intervention
    - Transplant prioritization

## Key Clinical Insights
Important predictors of mortality:
- Bilirubin → strongest risk indicator
- Ascites → ~3× higher mortality risk
- Prothrombin time → coagulation dysfunction
- Albumin → protective factor
These findings align with clinical knowledge and enhance model interpretability.

## Model Interpretability
- LASSO performs automatic feature selection
- RSF provides variable importance rankings
- Cox model offers clear hazard ratio interpretation
Balances Accuracy + clinical explanability

## Applications
- Clinical risk prediction
- Patient monitoring optimization
- Transplant timing decisions
- Healthcare resource allocation


