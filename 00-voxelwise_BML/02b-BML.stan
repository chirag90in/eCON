// The input data is a vector 'Y' of length 'N'.
data {
  int<lower=1> N;			// number of data points
  vector[N] Y;				// response variable
  vector[N] bpd;			// covariate
  int<lower=1> N_SUB;			// number of subjects
  int<lower=1> N_ROI;			// number of ROIs
  int<lower=1> N_VOX;			// number of VOX
  int<lower=1, upper=N_SUB> sid[N];	// subject id
  int<lower=1, upper=N_ROI> rid[N];	// ROI id
  int<lower=1, upper=N_VOX> vid[N];	// VOX id
}

// The parameters accepted by the model.
parameters {
  real POP;				  // fixed intercept
  real BPD;				  // fixed bpd slope
  vector[N_SUB] SUB;			  // subject varying intercept
  vector[N_ROI] ROI;			  // ROI varying intercept
  vector[N_VOX] VOX;			  // VOX varying intercept
  real<lower=0> sigma_e;		  // error sd
  real<lower=0> sigma_SUB;		  // subject sd
  real<lower=0> sigma_ROI;		  // ROI sd
  real<lower=0> sigma_VOX;		  // VOX sd
  real<lower=0> nu;			  // ddof
  real<lower=0> nu_SUB;			  // subject ddof
  real<lower=0> nu_ROI;			  // ROI ddof
  real<lower=0> nu_VOX; 		  // VOX ddof
}

// The model to be estimated.
model {
  real mu;
  // hyperpriors for standard deviation
  sigma_e ~ student_t(3, 0, 10);
  sigma_SUB ~ student_t(3, 0, 10);
  sigma_ROI ~ student_t(3, 0, 10);
  sigma_VOX ~ student_t(3, 0, 10);

  // hyperpriors for ddof
  nu ~ gamma(3.325,0.1);
  nu_SUB ~ gamma(3.325,0.1);
  nu_ROI ~ gamma(3.325,0.1);
  nu_VOX ~ gamma(3.325,0.1);
  
  // priors: assigning the same prior to all elements in case of vector parameters
  POP ~ student_t(3, 0, 10);		 // population intercept
  BPD ~ student_t(3, 0, 10);		 // fixed effect slope for bpd
  SUB ~ student_t(nu_SUB, 0, sigma_SUB); // subject intercepts
  ROI ~ student_t(nu_ROI, 0, sigma_ROI); // ROI intercepts
  VOX ~ student_t(nu_VOX, 0, sigma_VOX); // voxel intercepts
    
  // likelihood
  Y ~ student_t(nu, POP + BPD*bpd + SUB[sid] + ROI[rid] + VOX[vid], sigma_e);
}
