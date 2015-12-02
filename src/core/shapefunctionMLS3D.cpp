//-- Licencia --

#include "shapefunctionMLS3D.h"
#include "point.h"
#include "node.h"
#include "LMX/lmx.h"

namespace mknix {

ShapeFunctionMLS::ShapeFunctionMLS()
{
}


ShapeFunctionMLS::ShapeFunctionMLS( int nn_in,
                                    int mm_in,
                                    int weightType_in,
                                    double& alpha_c_in,
                                    double& d_c_in,
                                    Point* gp_in
                                  )
    : ShapeFunction(gp_in)
    , nn(nn_in)
    , mm(mm_in)
    , weightType(weightType_in)
{
    s_max = alpha_c_in * d_c_in;
//  cout << "S_max = " << s_max << endl;
    nn = this->gp->supportNodesSize;
    w.resize(6,nn);
    this->phi.resize(6, nn);
    if(mm == 1) {
        m = dim+1;
    }
    else std::cout << "WARNING: order not implemented." << std::endl;
    for(int i=0; i<m; ++i) {
        A.push_back( lmx::DenseMatrix<double>(m, m) );
        B.push_back( std::vector< lmx::Vector<double> >() );
        for (int j=0; j<nn; ++j) B[i].push_back( lmx::Vector<double>(m) );
    }

    P.resize(nn,m);
    p.resize(m);

//   std::ofstream sout("shapefunction.out");
//   sout << "dim = " << dim << endl;
//   sout << "nn (number of radial basis fun) = " << nn << endl;
//   sout << "mm (number of polynomial basis fun) = " << mm << endl;
//   sout << "rbfType = " << rbfType << endl;
//   sout << "alpha_c = " << alpha_c << endl;
//   sout << "d_c = " << d_c << endl;
//   sout << "q = " << q << endl;
//   sout << "gp = " << gp->i << endl;

}


ShapeFunctionMLS::~ShapeFunctionMLS()
{
}

void ShapeFunctionMLS::calc()
{
    computeWeights();
    computeMomentMatrix();
    // Solve linear eq for shape function:
    computePhi(gp->X, gp->Y);
//   computePhi(0.2, 0.4);
}

void ShapeFunctionMLS::computeWeights()
{
//  double radius;
    double rr2, xn, yn, zn(0), xp, yp, zp(0), qq;
    //For case 1:
    double s_i, ss_i;


//  if (dim != 2)
//    std::cout << "Warning, space dimension is not implemented." << std::endl;
//  else{
    // Switch type of weight function:
    switch (weightType) {
    case 0:
        // weightType = 0 -> EXP
        // Compute w(i) = exp( -( (x_i - x_j)^2 + (y_i - y_j)^2 )
        //                / (s_max*alpha)^2 )
        //              = exp( -rr2 * qq )
        qq = 1. / std::pow(s_max*0.4, 2);
        xp = gp->X;
        yp = gp->Y;
	if(dim==3) zp = gp->Z;
        for (int j=0; j<nn; ++j) {
            xn = gp->supportNodes[j]->getX();
            yn = gp->supportNodes[j]->getY();
	    if(dim==3) zn = gp->supportNodes[j]->getZ();
            rr2 = std::pow( xp - xn, 2) +
                  std::pow( yp - yn, 2);
	    if(dim==3) rr2 += std::pow( zp - zn, 2);
            w(0,j) = std::exp( -qq * rr2);
//          cout << s_max << ": " << rr2 << " ";
//          cout << w(0,j)<< endl;
            w(1,j) = -2 * qq * std::exp( -qq * rr2) * ( xp - xn );
            w(2,j) = -2 * qq * std::exp( -qq * rr2) * ( yp - yn );
	    if(dim==3)
	      w(3,j) = -2 * qq * std::exp( -qq * rr2) * ( zp - zn );
	    // for plates and shells:
//             w(3,j) = -2 * qq * std::exp( -qq * rr2)
//                      + 4 * qq * qq * std::pow( xp - xn, 2 ) * std::exp( -qq * rr2);
//             w(4,j) = 4 * qq * qq * std::exp( -qq * rr2)
//                      * ( xp - xn ) * ( yp - yn );
//             w(5,j) = -2 * qq * std::exp( -qq * rr2)
//                      + 4 * qq * qq * std::pow( yp - yn, 2 ) * std::exp( -qq * rr2);
        }
        break;

    case 1:
        // weightType = 0 -> Cubic Spline
        // Compute w(i) =
        xp = gp->X;
        yp = gp->Y;
	if(dim==3) zp = gp->Z;
        for (int j=0; j<nn; ++j) {
            xn = gp->supportNodes[j]->getX();
            yn = gp->supportNodes[j]->getY();
	    if(dim==3) zn = gp->supportNodes[j]->getZ();
            if(dim==2) s_i = std::sqrt(std::pow( xn - xp, 2) +
                                       std::pow( yn - yp, 2) );
	    if(dim==3) s_i = std::sqrt(std::pow( xn - xp, 2) +
                                       std::pow( yn - yp, 2) +
				       std::pow( zn - zp, 2) );
            ss_i = s_i/s_max;

            if(ss_i < .5) {
                w(0,j) = 2./3. - 4.*pow(ss_i,2) + 4.*pow(ss_i,3);
                w(1,j) = (8. - 12.*ss_i)*(xn-xp)/(s_max*s_max);
                w(2,j) = (8. - 12.*ss_i)*(yn-yp)/(s_max*s_max);
		if(dim==3) w(3,j) = (8. - 12.*ss_i)*(zn-zp)/(s_max*s_max);
            }
            else if (ss_i <= 1.) {
                w(0,j) = 4./3. - 4*ss_i + 4.*pow(ss_i,2) - (4./3.)*pow(ss_i,3);
                w(1,j) = ( 4./ss_i - 8. + 4.*ss_i)*(xn-xp)/(s_max*s_max);
                w(2,j) = ( 4./ss_i - 8. + 4.*ss_i)*(yn-yp)/(s_max*s_max);
		if(dim==3) w(3,j) = ( 4./ss_i - 8. + 4.*ss_i)*(zn-zp)/(s_max*s_max);
            }
//          cout << s_max << ": " << rr2 << " ";

//          cout << "(" << xn << ","<< yn <<") : "
//               << w(0,j) <<", "
//               << w(1,j) <<", "
//               << w(2,j) <<", "
//               << endl;
        }
        break;

    }
//  }

    // Polynomials of weight points:
    //P_{ij} = P_i(X_j - X)
    if(mm == 1) {
        for (int i=0; i<nn; ++i) {
            P(i,0) = 1.;
            P.writeElement( gp->supportNodes[i]->getX(), i, 1);
            P.writeElement( gp->supportNodes[i]->getY(), i, 2);
            if( dim == 3 ) {
                P.writeElement( gp->supportNodes[i]->getZ(), i, 3);
            }
        }
    }
    else
        std::cout
                << "WARNING: order not implemented in polynomials of weight functions."
                << std::endl;
}

void ShapeFunctionMLS::computeMomentMatrix()
{
    // Example with mm=1 (first order):
    // [A] = P^T W P
    // Dimensions in 2D...
    // [A]_{3x3} = [P^T]_{3 x nn} W_{nn x nn} P_{nn x 3}
    // index mult:
    // A_{ij} = w_k*p_i(x_k)*p_j(x_k)

    if(mm == 1) {
        // A = sum_i w(i)*ppt(i)
        //  with:
        //  ppt(i) = p*p^T = [[1 xi yi] [xi xi^2 xi*yi] [yi xi*yi yi^2]]
        for (int i=0; i<m; ++i) {
            for (int j=0; j<m; ++j) {
                for (int k=0; k<nn; ++k) {
                    A[0].addElement( w.readElement(0,k) *
                                     P.readElement(k,i) *
                                     P.readElement(k,j)
                                     , i, j
                                   );
//         cout << w.readElement(0,k) << " * "
//              << P.readElement(i,k) << " * "
//              << P.readElement(j,k) << " = "
//              << A[0](i,j) << endl;
                }
            }
        }
        //DERIVATIVES FOR X and Y... only for mm=1 (vamos, chapucilla...)
        // A,x = sum_i ( w(i),x * ppt(i) )
        for (int i=0; i<m; ++i) {
            for (int j=0; j<m; ++j) {
                for (int k=0; k<nn; ++k) {
                    A[1].addElement( w.readElement(1,k) * // derivada x de weight(k)
                                     P.readElement(k,i) *
                                     P.readElement(k,j)
                                     , i, j
                                   );
                }
            }
        }
        // A,y = sum_i ( w(i),y * ppt(i) )
        for (int i=0; i<m; ++i) {
            for (int j=0; j<m; ++j) {
                for (int k=0; k<nn; ++k) {
                    A[2].addElement( w.readElement(2,k) * // derivada y de weight(k)
                                     P.readElement(k,i) *
                                     P.readElement(k,j)
                                     , i, j
                                   );
                }
            }
        }
        if( dim == 3 ) {
        // A,z = sum_i ( w(i),z * ppt(i) )
	  for (int i=0; i<m; ++i) {
	      for (int j=0; j<m; ++j) {
		  for (int k=0; k<nn; ++k) {
		      A[3].addElement( w.readElement(3,k) * // derivada z de weight(k)
				      P.readElement(k,i) *
				      P.readElement(k,j)
				      , i, j
				    );
		  }
	      }
	  }
	}
        // Now we compute the RHS: B_I
        // => B_I(i) = P^T * W = w_I * P_Ii
        // And the derivatives of the RHS: B_I,x
        // => B_I(i) = P^T * W,x = w_I,x * P_Ii
        for (int j=0; j<nn; ++j) {
            for (int i=0; i<m; ++i) {
                for (int r=0; r<m; ++r) {
                    B[r][j].writeElement( w.readElement(r,j) * P.readElement(j,i), i );
                }
            }
        }
    }
    else std::cout << "WARNING: order not implemented." << std::endl;
}

void ShapeFunctionMLS::computePhi(double xp, double yp, double zp)
{
    std::vector< lmx::Vector<double> > alpha; /**< auxiliary variables */
//  lmx::Vector<double> alpha_1(m); /**< auxiliary variables, x derivative */
//  lmx::Vector<double> alpha_2(m); /**< auxiliary variables, y derivative */


    // Solving linear system A * alpha_0 = p(gp) = rhs;
    // PHI and derivatives (6 systems are solved):
    //////////////////////////////////////
    // CHANGED FOR SOLVING ONLY 3 SYSTEMS: PHI AND FIRST PARTIAL
    // DERIVATIVES (dX AND dY)
    //////////////////////////////////////
    // PHI:
    {
        alpha.push_back(lmx::Vector<double>(m));
        lmx::Vector<double> rhs(m);
        rhs(0) = 1.0; // base is shifted to GP, so rest of elements are zero.
        rhs(1) = gp->X; // base is shifted to GP, so rest of elements are zero.
        rhs(2) = gp->Y; // base is shifted to GP, so rest of elements are zero.
        if(dim == 3)
	  rhs(3) = gp->Z; 
        lmx::LinearSystem<double> lsys(A[0], alpha[0], rhs);
        lsys.solveYourself();
        for (int j = 0; j < nn; ++j) {
            // phi_I = alpha * B_I
            this->phi(0,j) = alpha[0]*B[0][j];
        }
    }

    if(mm == 1) {
        alpha.push_back(lmx::Vector<double>(m)); //x derivative
        alpha.push_back(lmx::Vector<double>(m)); //y derivative
        if(dim == 3)
	  alpha.push_back(lmx::Vector<double>(m)); //z derivative
        // first the x derivative:
        {
            // A*alpha,x = p,x  - A,x*alpha
            lmx::Vector<double> rhs(m);
            rhs -= A[1] * alpha[0]; // = -A,x*alpha
            rhs.addElement( 1. , 1 ); // += p,x = [0 1 0]
            lmx::LinearSystem<double> lsys(A[0], alpha[1], rhs);
            lsys.solveYourself();
            // phi,x
            for (int j = 0; j < nn; ++j) {
                this->phi(1,j)
//        = alpha[1] * B[0][j] + alpha[0] * B[1][j];
                += (alpha[0](0) + alpha[0](1)*P(j,1) + alpha[0](2)*P(j,2))*w(1,j)
                   + (alpha[1](0) + alpha[1](1)*P(j,1) + alpha[1](2)*P(j,2))*w(0,j);
//        for(int k=0;k<m;++k){
//          cout << "phi(1,"<<j<<") +="
//          << alpha[1](k) << "*"<< B[0][j](k) << "+"
//          << alpha[0](k) << "*"<< B[1][j](k) << endl;
//        }
//        cout << "phi(1,"<<j<<") = " << phi(1,j) << endl;
            }
        }
        // the y derivative:
        {
            // A*alpha,y = p,y  - A,y*alpha
            lmx::Vector<double> rhs(m);
            rhs -= A[2] * alpha[0]; // = -A,y*alpha
            rhs.addElement( 1. , 2 ); // += p,y = [0 0 1]
            lmx::LinearSystem<double> lsys(A[0], alpha[2], rhs);
            lsys.solveYourself();
            // phi,y
            for (int j = 0; j < nn; ++j) {
                this->phi(2,j)
                += (alpha[0](0) + alpha[0](1)*P(j,1) + alpha[0](2)*P(j,2))*w(2,j)
                   + (alpha[2](0) + alpha[2](1)*P(j,1) + alpha[2](2)*P(j,2))*w(0,j);
            }
        }
        // the z derivative:
        if(dim == 3){
            // A*alpha,z = p,z  - A,z*alpha
            lmx::Vector<double> rhs(m);
            rhs -= A[3] * alpha[0]; // = -A,z*alpha
            rhs.addElement( 1. , 3 ); // += p,z = [0 0 0 1]
            lmx::LinearSystem<double> lsys(A[0], alpha[3], rhs);
            lsys.solveYourself();
            // phi,z
            for (int j = 0; j < nn; ++j) {
                this->phi(3,j)
                += (alpha[0](0) + alpha[0](1)*P(j,1) + alpha[0](2)*P(j,2))*w(3,j)
                   + (alpha[3](0) + alpha[3](1)*P(j,1) + alpha[3](2)*P(j,2))*w(0,j);
            }
        }
    }

    //Clean possible round off errors..
    phi.clean(1E-28);
    double sumphi = 0;
    for(int i=0; i<nn ; ++i) {
        sumphi += this->phi(0,i);
    }
//  cout << "nPs: " << nn << ", SUM_PHI = " << sumphi << endl;
//   cout << "PHI = " << phi << endl;

//  this->outputValues();
//  this->gnuplotOut();
}

} //Namespace mknix
