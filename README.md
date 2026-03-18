# Survival-Modeling-of-Cirrhosis-Patients-Using-Statistical-Learning-Methods

## Overview
### 
This project develops and evaluates machine learning–based survival models to predict mortality risk in patients with liver cirrhosis. The goal is to support clinical decision-making, including patient monitoring, risk stratification, and transplant timing.
Using real-world clinical data from a Mayo Clinic study, we compare traditional statistical approaches with modern machine learning techniques to improve predictive performance and interpretability.

## Dataset
###
Source: Mayo Clinic cirrhosis study (Primary Biliary Cirrhosis)
Dataset link (Kaggle): (https://www.kaggle.com/datasets/joebeachcapital/cirrhosis-patient-survival-prediction)
Sample size: 418 patients
Features: 17 clinical variables (demographics, symptoms, lab biomarkers)
Outcome:
Death
Censored (alive)
Liver transplant (treated as competing risk)

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

- 
