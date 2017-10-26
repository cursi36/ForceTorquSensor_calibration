
# INSTRUCTIONS FOR USING CODES FOR LIVE CALIBRATION

- Run on terminal

 ~/ecat_dev/ec_master_test.git/build$ ati_ft6_calib ../configs/nrt_config.yaml

(the file reads from sensor and saves sens.txt on runtime. Choose proper folder for saving sens.txt).

- Open /home/user/ecat_dev/ec_master_test.git/examples/ ati_ft6_calib/ft_ati_test.cpp to change file to save.

- Run it without loading the sensor so as to get load-free sensor's data.

- Open Matlab and run acquire_data.m
- Choose Calibrate. This runs LiveCalibrate.m for live calibration.
- Change LiveCalibrate.m to change constraints or minimization methods.

