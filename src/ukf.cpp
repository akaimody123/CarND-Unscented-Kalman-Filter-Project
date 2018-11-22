#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);


  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5; //original 30

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5; //original 30
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  /////////////////////////////////////////////////////////////////
 
  //newly added
  
  n_x_=5;

  ///* Augmented state dimension
  n_aug_ = n_x_ + 2;

  //lambda for augmented states
  lambda_ = 3 - n_aug_;

  ///* Sigma points dimension
  n_sig_ = 2 * n_aug_ + 1;

  // Initialize weights.
  weights_ = VectorXd(n_sig_);
  weights_.fill(0.5 / (n_aug_ + lambda_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  // Initialize measurement noice covarieance matrix
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0,std_radrd_*std_radrd_;

  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << std_laspx_*std_laspx_,0,
              0,std_laspy_*std_laspy_;
  
////////////////////////////////////////////////////////////
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if(!is_initialized_)
  {
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) //radar measurement
  { double rho = meas_package.raw_measurements_(0);
    double phi = meas_package.raw_measurements_(1);
    double rho_dot = meas_package.raw_measurements_(2);
    x_ << rho*cos(phi),rho*sin(phi),rho_dot, 0, 0;
    std::cout << "radar initialized ok"  << std::endl;
	
  }
  else
  { x_ << meas_package.raw_measurements_(0),meas_package.raw_measurements_(1),0, 0, 0; //lidar measurement
    std::cout << "lidar initialized ok"  << std::endl;
  }

  // Saving first timestamp in seconds
  time_us_ = meas_package.timestamp_ ;
  // done initializing, no need to predict or update
  is_initialized_ = true;

  return;
  }

  // Calculate dt
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  // Prediction step
  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  }
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  }



}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  //initialze augmented states and covaranice matrix
  x_aug = VectorXd(n_aug_); //augmented state
  x_aug.head(5)=x_;
  x_aug(5)=0; //process noise has 0 mean
  x_aug(6)=0;

  P_aug=MatrixXd(n_aug_,n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_)=P_;
  P_aug(5,5)=std_a_*std_a_;
  P_aug(6,6)=std_yawdd_*std_yawdd_;

  std::cout << "x_aug,P_aug ini ok"  << std::endl;

  //generate sigma points
  x_sig = MatrixXd(n_aug_,n_sig_);

  A=MatrixXd(n_aug_,n_aug_);

  A=P_aug.llt().matrixL();

  std::cout << "A ok"  << std::endl;

  x_sig.col(0)=x_aug;

  std::cout << "x_sig.col(0) ok"  << std::endl;

  for(int i=0; i<n_aug_;i++)
  	{
  	 	x_sig.col(i+1)=x_aug+sqrt(lambda_+n_aug_)*A.col(i);
  	
  	 	x_sig.col(i+1+n_aug_)=x_aug-sqrt(lambda_+n_aug_)*A.col(i);
  	}
  std::cout << "x_sig"  << x_sig << std::endl;

  //predict sigma points--extract info and calculate state&covarance matrix in CTRV model
  //predict mean

  XsigPred=MatrixXd(n_x_,n_sig_);
  XsigPred.fill(0.0);

  for (int i=0;i<n_sig_;i++)
  {
  	//extract state info
  	float px=x_sig(0,i);
  	float py=x_sig(1,i);
  	float v=x_sig(2,i);
  	float yaw=x_sig(3,i);
  	float yawd=x_sig(4,i);
  	float nu_a=x_sig(5,i);
  	float nu_yawdd=x_sig(6,i);

  	//calculate mean using CTRV model
  	xsig_pred=VectorXd(n_x_);
  	xdiff=VectorXd(n_x_);
  	nu_diff=VectorXd(n_x_);
  	if (yawd > 0.001 || yawd<-0.001)
  	{
  		xdiff<< v/yawd*(sin(yaw+yawd*delta_t)-sin(yaw)),v/yawd*(-cos(yaw+yawd*delta_t)+cos(yaw)),0,yawd*delta_t,0;
  		nu_diff<<0.5*delta_t*delta_t*cos(yaw)*nu_a,0.5*delta_t*delta_t*sin(yaw)*nu_a,delta_t*nu_a,0.5*delta_t*delta_t*nu_yawdd,delta_t*nu_yawdd;
  		//xsig_pred=x_sig(0,i,n_x_,1)+xdiff+nu_diff;//x_+xdiff+nu_diff;

  	}
  	else
  	{
  		xdiff<<v*cos(yaw)*delta_t,v*sin(yaw)*delta_t,0,yawd*delta_t,0;
  		nu_diff<<0.5*delta_t*delta_t*cos(yaw)*nu_a,0.5*delta_t*delta_t*sin(yaw)*nu_a,delta_t*nu_a,0.5*delta_t*delta_t*nu_yawdd,delta_t*nu_yawdd;
  		//xsig_pred=x_sig(0,i,n_x_,1)+xdiff+nu_diff;//x_+xdiff+nu_diff;

  	}
  	
  	//XsigPred.col(i)=xsig_pred;
  XsigPred(0,i)=px+xdiff(0)+nu_diff(0);
	XsigPred(1,i)=py+xdiff(1)+nu_diff(1);
	XsigPred(2,i)=v+xdiff(2)+nu_diff(2);
	XsigPred(3,i)=yaw+xdiff(3)+nu_diff(3);
	XsigPred(4,i)=yawd+xdiff(4)+nu_diff(4);

  }
 // std::cout << "XsigPred"  << XsigPred << std::endl;


  	//predict mean
    XsigPred_mean=VectorXd(n_x_);
    XsigPred_mean.fill(0.0);
 

  	for(int i=0;i<n_sig_;i++)
  	{
  		
  		XsigPred_mean+=weights_(i)*XsigPred.col(i);
  	} 


 


    //normalize angle
    while (XsigPred_mean(3)>M_PI)
    {
    	XsigPred_mean(3)-=2*M_PI;
    }
    while (XsigPred_mean(3)<-M_PI)
    {
    	XsigPred_mean(3)+=2*M_PI;
    }

  	

  	//predict covariance
    XsigPred_cov=MatrixXd(n_x_,n_x_);

    XsigPred_cov.fill(0.0);

    XsigPred_diff=VectorXd(n_x_);
    

   	for(int i=0;i<n_sig_;i++)
  	{
		//normalize angle

		XsigPred_diff=XsigPred.col(i)-XsigPred_mean;
    		while (XsigPred_diff(3)>M_PI)
    		{XsigPred_diff(3)-=2*M_PI;}
    		while (XsigPred_diff(3)<-M_PI)
    		{XsigPred_diff(3)+=2*M_PI;}
  		
  		XsigPred_cov+=weights_(i)*XsigPred_diff*XsigPred_diff.transpose();
                
  	} 

  	//assign predicted state and covarance matrix
  	x_=XsigPred_mean;

  	P_=XsigPred_cov;

   // std::cout << "XsigPred_mean " << XsigPred_mean << std::endl;
   // std::cout << "XsigPred_cov " << XsigPred_cov << std::endl;

  







 }



/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
	//predict measurement

  	//calculate new sigma points in measurement model
  	Zsig_pred=MatrixXd(2,n_sig_);
  	Zsig_pred.fill(0.0);

  	Zsig_pred_mean=VectorXd(2);
  	Zsig_pred_mean.fill(0.0);

  	Zsig_pred_cov=MatrixXd(2,2);
  	Zsig_pred_cov.fill(0.0);
  	if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  	{
  		for(int i=0; i<n_sig_;i++)
  		{
  			Zsig_pred.col(i)<<XsigPred(0,i),XsigPred(1,i);

  		}

  		//calculate predicted measurement mean
  		for(int i=0;i<n_sig_;i++)
  		{
  			Zsig_pred_mean+=weights_(i)*Zsig_pred.col(i);
  		}
  		//calculate predicted measurement covariance
  		for(int i=0;i<n_sig_;i++)
  		{
  			Zsig_pred_cov+=weights_(i)*(Zsig_pred.col(i)-Zsig_pred_mean)*(Zsig_pred.col(i)-Zsig_pred_mean).transpose();
  		}

  		Zsig_pred_cov+=R_lidar_;

      std::cout << "Zsig_pred_mean_lidar " << Zsig_pred_mean << std::endl;
      std::cout << "Zsig_pred_cov_lidar " << Zsig_pred_cov << std::endl;

  		//UKF update

  		//correlation between state space and measurement space
  		T=MatrixXd(n_x_,2);
  		T.fill(0.0);
  		for(int i=0;i<n_sig_;i++)
  		{
  			T+=weights_(i)*(XsigPred.col(i)-x_)*(Zsig_pred.col(i)-Zsig_pred_mean).transpose();
  		}

    //  std::cout << "T_lidar " << T << std::endl;

  		//kalman gain
  		K=MatrixXd(n_x_,2);
  		K=T*Zsig_pred_cov.inverse();

    //  std::cout << "K_lidar " <<K << std::endl;

  		//state update
  		z=VectorXd(2);
  		z=meas_package.raw_measurements_;

                ZsigPred_diff=z-Zsig_pred_mean;

  		x_=x_+K*ZsigPred_diff;

    //  std::cout << "x_lidar " << x_ << std::endl;

  		//covariance update
  		P_=P_-K*Zsig_pred_cov*K.transpose();

    //  std::cout << "P_lidar " << P_ << std::endl;

                 //calculate NIS
      NIS_lidar_ = ZsigPred_diff.transpose() * Zsig_pred_cov.inverse() * ZsigPred_diff;
      std::cout << "NIS_lidar " << NIS_lidar_ << std::endl;
}


}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */


	//predict measurement

  	//calculate new sigma points in measurement model
  	Zsig_pred=MatrixXd(3,n_sig_);
  	Zsig_pred.fill(0.0);

  	Zsig_pred_mean=VectorXd(3);
  	Zsig_pred_mean.fill(0.0);

  	Zsig_pred_cov=MatrixXd(3,3);
  	Zsig_pred_cov.fill(0.0);

	ZsigPred_diff=VectorXd(3);

  	if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  	{
  		for(int i=0; i<n_sig_;i++)
  		{
  			Zsig_pred.col(i)<<sqrt(XsigPred(0,i)*XsigPred(0,i)+XsigPred(1,i)*XsigPred(1,i)), atan2(XsigPred(1,i),XsigPred(0,i)),(XsigPred(0,i)*XsigPred(2,i)*cos(XsigPred(3,i))+XsigPred(1,i)*XsigPred(2,i)*sin(XsigPred(3,i)))/sqrt(XsigPred(0,i)*XsigPred(0,i)+XsigPred(1,i)*XsigPred(1,i));

  		}

  		//calculate predicted measurement mean
  		for(int i=0;i<n_sig_;i++)
  		{
  			Zsig_pred_mean+=weights_(i)*Zsig_pred.col(i);
  		}

		//normalize angle
    		while (Zsig_pred_mean(1)>M_PI)
    		{Zsig_pred_mean(1)-=2*M_PI;}
    		while (Zsig_pred_mean(1)<-M_PI)
    		{Zsig_pred_mean(1)+=2*M_PI;}

  		//calculate predicted measurement covariance
  		for(int i=0;i<n_sig_;i++)
  		{

			ZsigPred_diff=Zsig_pred.col(i)-Zsig_pred_mean;
	    		while (ZsigPred_diff(1)>M_PI)
	    		{ZsigPred_diff(1)-=2*M_PI;}
	    		while (ZsigPred_diff(1)<-M_PI)
	    		{ZsigPred_diff(1)+=2*M_PI;}
  			Zsig_pred_cov+=weights_(i)*ZsigPred_diff*ZsigPred_diff.transpose();
  		}

  		Zsig_pred_cov+=R_radar_;

      // std::cout << "Zsig_pred_mean_radar " << Zsig_pred_mean << std::endl;
      // std::cout << "Zsig_pred_cov_radar " << Zsig_pred_cov << std::endl;


  		//UKF update

  		//correlation between state space and measurement space
  		T=MatrixXd(n_x_,3);
  		T.fill(0.0);
  		for(int i=0;i<n_sig_;i++)
  		{
			ZsigPred_diff=Zsig_pred.col(i)-Zsig_pred_mean;
	    		while (ZsigPred_diff(1)>M_PI)
	    		{ZsigPred_diff(1)-=2*M_PI;}
	    		while (ZsigPred_diff(1)<-M_PI)
	    		{ZsigPred_diff(1)+=2*M_PI;}

			XsigPred_diff=XsigPred.col(i)-x_;
	    		while (XsigPred_diff(3)>M_PI)
	    		{XsigPred_diff(3)-=2*M_PI;}
	    		while (XsigPred_diff(3)<-M_PI)
	    		{XsigPred_diff(3)+=2*M_PI;}

			
  			T+=weights_(i)*XsigPred_diff*ZsigPred_diff.transpose();
  		}

      //std::cout << "T_radar " << T << std::endl;


  		//kalman gain
  		K=MatrixXd(n_x_,3);
  		K=T*Zsig_pred_cov.inverse();

      //std::cout << "K_radar " << K << std::endl;


  		//state update
  		z=VectorXd(3);
  		z=meas_package.raw_measurements_;

		ZsigPred_diff=z-Zsig_pred_mean;
    		while (ZsigPred_diff(1)>M_PI)
    		{ZsigPred_diff(1)-=2*M_PI;}
    		while (ZsigPred_diff(1)<-M_PI)
    		{ZsigPred_diff(1)+=2*M_PI;}

  		x_=x_+K*ZsigPred_diff;

  		//covariance update
  		P_=P_-K*Zsig_pred_cov*K.transpose();

     // std::cout << "x_radar " << x_ << std::endl;
     // std::cout << "P_radar " << P_ << std::endl;

      //calculate NIS
      NIS_radar_ = ZsigPred_diff.transpose() * Zsig_pred_cov.inverse() * ZsigPred_diff;
      std::cout << "NIS_radar " << NIS_radar_ << std::endl;

 


  	}


}

