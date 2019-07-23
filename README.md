start with cali_tutorial.R

Logistic Regression and Support Vector Machine are trained in cali_tutorial.R

AUC is measured with built in package

calibration is measured with:<br/>
spiegelhalter Z test: Spiegelhalter_z.R<br/>
Homser Lemeshow C or H test: hosmer_lemeshow.R<br/>
Expected calibration error and Maximum calibration error: mce_ece.R<br/>
reliability diagram: reliability_diagram.R<br/>

recalibration is done with:<br/>
Platt scaling: lr_platt_recal.R and svm_platt_recal.R<br/>
isotonic regression: lr_iso_recal.R and svm_iso_recal.R<br/>
