This is a code repository for the following book chapter:

Chang, W, Haran, M. (2024+), Computer Model Emulation and Calibration for Discrete Data. In Bingham, D., Haran, M., Oakley, J., Sanso, B. (Ed.), Handbook of Statistical Methods for Computer Modeling
Please contact the first author at wonchang@snu.ac.kr for any questions.


SIRbayesian_GLM.R: Code for calibrating SIR model example
VSIRbayesian_GLM.R: Code for calibrating SIR model with a vaccination effect example

The following code files are for calibrating the ice model calibration using semi-continuous spatial data example:
1_read_data.R: Preparing data objects for emulation
2_emulator.R: GP-based emulator for semi-continuous spatial data
2-1_emulator_check.R: Cross-validation to check the emulation performance
3_calibration.R: Bayesian calibration using the emulator created by '2_emulator.R'


The methodology is originally developed by

Chang, W., Konomi, B. A., Karagiannis, G., Guan, Y., Haran, M. (2022) Ice model calibration using semi-continuous spatial data, the Annals of Applied Statistics, 16 (3), 1937-1961
