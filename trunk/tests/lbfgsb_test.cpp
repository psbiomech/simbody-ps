#include "Simmath.h"
#include "SimTKcommon.h"
#include "SimTKcommon/internal/common.h"
#include "simmatrix/internal/BigMatrix.h"
#include "Optimizer.h"

#include <iostream>
using std::cout;
using std::endl;


const static int NUMBER_OF_PARAMETERS = 25;

class ProblemSystem : public SimTK::OptimizerSystem {

   public:

   ProblemSystem( const int numParameters ) : SimTK::OptimizerSystem( numParameters ) {}

   int objectiveFunc(   SimTK::Vector &coefficients, bool new_coefficients,  double *f  ) const  {
      double *x;
      int i;

      x = &coefficients[0];

//printf("objectiveFunction x = ",x[0],x[1],x[2]);
      *f = .25 *(x[0]-1.0)*(x[0]-1.0);
//   printf(" %f",x[0]);
      for(i=1;i<numParameters;i++) {
         *f = *f + pow(x[i]-x[i-1]*x[i-1], 2.0);
//   printf(" %f",x[i]);
      }

//   printf(" \n");
      *f = 4.0* *f;
      return( 0 ); 
   }

   int gradientFunc( SimTK::Vector &coefficients, bool new_coefficients,  SimTK::Vector &gradient ) const {
      double *x,t1,t2;
      int i;

      x = &coefficients[0]; 

      t1 = x[1]-(x[0]*x[0]);
      gradient[0] = 2.0*(x[0]-1.0)-16.0*x[0]*t1;
      for(i=1;i<numParameters-1;i++) {
         t2=t1;
         t1=x[i+1]-(x[i]*x[i]);
         gradient[i]=8.0*t2-16.0*x[i]*t1;
      }
      gradient[numParameters-1]=8.0*t1;
// printf("objectiveGradient x = %f %f %f  g = %f \n",x[0],x[1],x[2],gradient[0]);

    return(0);

   }

};

/* adapted from driver1.f of Lbfgsb.2.1.tar.gz  */
main() {

    double params[10],f;
    int i;
    int n = NUMBER_OF_PARAMETERS;

    SimTK::Vector results(NUMBER_OF_PARAMETERS);
    SimTK::Vector lower_bounds(NUMBER_OF_PARAMETERS);
    SimTK::Vector upper_bounds(NUMBER_OF_PARAMETERS);

    ProblemSystem sys(NUMBER_OF_PARAMETERS);

    cout << "LBFGSB driver1 test " << endl;

    /* set initial conditions */
    for(i=0;i<n;i++) {
       results[i] = 3.0;
    }

    /* set bounds */
    for(i=0;i<n;i=i+2) {   // even numbered 
       lower_bounds[i] = 1.0;
       upper_bounds[i] = 100.0;
    }
    for(i=1;i<n;i=i+2) { // odd numbered
       lower_bounds[i] = -100.0;
       upper_bounds[i] = 100.0;
    }

    sys.setParameterLimits( lower_bounds, upper_bounds );

    try {
    SimTK::Optimizer opt( sys ); 

    params[0] = 100;
    opt.setOptimizerParameters( MAX_FUNCTION_EVALUATIONS, params );

    params[0] = .0001;
    opt.setOptimizerParameters( GRADIENT_CONVERGENCE_TOLERANCE, params );

    params[0] = 1.0;
    opt.setOptimizerParameters( DEFAULT_STEP_LENGTH, params );

    params[0] = 0.9;
    opt.setOptimizerParameters( LINE_SEARCH_ACCURACY, params );


    f = opt.optimize( results );

    }

    catch (SimTK::Exception::Base exp) {
        cout << "Caught exception :" << exp.getMessage() << endl;
    }

    printf("f = %f params = ",f);
    for( i=0; i<NUMBER_OF_PARAMETERS; i++ ) {
       printf(" %f",results[i]); 
    }
    printf("\n");

}
