#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

#define LASER_COV         0.0225f
#define RADAR_RO_COV      0.09f
#define RADAR_PHI_COV     0.0009f
#define RADAR_RO_DOT_COV  0.09f

#define NOISE_AX       9
#define NOISE_AY       9

#define P_POS_INIT     10
#define P_VEL_INIT     1000

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  /**
  TODO:
    * Finish initializing the FusionEKF.
  */

  R_laser_ << LASER_COV, 0.0,
              0.0, LASER_COV;

  R_radar_ << RADAR_RO_COV, 0.0, 0.0,
              0.0, RADAR_PHI_COV, 0.0,
              0.0, 0.0, RADAR_RO_DOT_COV;

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << P_POS_INIT, 0, 0, 0,
             0, P_POS_INIT, 0, 0,
             0, 0, P_VEL_INIT, 0,
             0, 0, 0, P_VEL_INIT;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    //cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

      //set the state with the initial location and zero velocity
      float ro = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      ekf_.x_ << ro*cos(phi), ro*sin(phi), 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */

      //set the state with the initial location and zero velocity
      float px = measurement_pack.raw_measurements_[0];
      float py = measurement_pack.raw_measurements_[1];
      ekf_.x_ << px, py, 0, 0;
    }

    // initialize time
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
   */
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;
  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  cout << "dt = " << dt << endl;
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  //set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*NOISE_AX, 0, dt_3/2*NOISE_AX, 0,
             0, dt_4/4*NOISE_AY, 0, dt_3/2*NOISE_AY,
             dt_3/2*NOISE_AX, 0, dt_2*NOISE_AX, 0,
             0, dt_3/2*NOISE_AY, 0, dt_2*NOISE_AY;


  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    cout << "radar update..." << endl;
    ekf_.R_ = R_radar_; 
    Tools tools;
    bool is_jac_valid = false;
    Hj_ = tools.CalculateJacobian(ekf_.x_, is_jac_valid);
    if (is_jac_valid) {
      ekf_.H_ = Hj_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    }
  } else {
    // Laser updates
    cout << "laser update..." << endl;
    ekf_.R_ = R_laser_; 
    ekf_.H_ = H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
