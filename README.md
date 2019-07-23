start with cali_tutorial.R

Logistic Regression and Support Vector Machine are trained in cali_tutorial.R

AUC is measured with built in package

calibration is measured with:
spiegelhalter Z test: Spiegelhalter_z.R
Homser Lemeshow C or H test: hosmer_lemeshow.R
Expected calibration error and Maximum calibration error: mce_ece.R
reliability diagram: reliability_diagram.R

recalibration is done with:
Platt scaling: lr_platt_recal.R and svm_platt_recal.R
isotonic regression: lr_iso_recal.R and svm_iso_recal.R
