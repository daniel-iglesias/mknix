//-- Licencia --

#include "shapefunctionRBF.h"
#include "point.h"
#include "node.h"
#include "LMX/lmx.h"

namespace mknix {

ShapeFunctionRBF::ShapeFunctionRBF()
{
}


ShapeFunctionRBF::ShapeFunctionRBF( int nn_in,
                                    int mm_in,
                                    int rbfType_in,
                                    double& alpha_c_in,
                                    double& d_c_in,
                                    double& q_in,
                                    Point* gp_in
                                  )
    : ShapeFunction(gp_in)
    , nn(nn_in)
    , mm(mm_in)
    , rbfType(rbfType_in)
    , alpha_c(alpha_c_in)
    , d_c(d_c_in)
    , q(q_in)
{
    g_o.resize(nn + mm, nn + mm);
    if (dim == 2) {
        this->phi.resize(6, nn + mm);
    } else if (dim == 3) {
        std::cout << "Warning, space dimension is partially implemented."
                  << std::endl;
        this->phi.resize(4, nn+mm);
    }

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


ShapeFunctionRBF::~ShapeFunctionRBF()
{
}

void ShapeFunctionRBF::calc()
{
    computeMomentMatrix();
    // Solve linear eq for shape function:
    computePhi(gp->getX(), gp->getY(), 0);
//     computePhi(gp->getX(), gp->getY(), gp->getZ());
//   computePhi(0.2, 0.4);
}

void ShapeFunctionRBF::computeMomentMatrix()
{
    double radius, qq;

    if (dim != 2) {
        std::cout << "Warning, space dimension is partially implemented."
                  << std::endl;
    }
    // Switch type of RBF:
    switch (rbfType) {
    case 0:
        // rbfType = 0 -> MQ
        // Compute R(i,j) = ((x_i - x_j)^2 + (y_i - y_j)^2 + (alpha_c*d_c)^2 )^q
        for (auto i = 0u; i < nn; ++i) {
            for (auto j = 0u; j < nn; ++j) {
                g_o(i, j) = std::pow(std::pow(gp->supportNodes[i]->getX() -
                                              gp->supportNodes[j]->getX(), 2)
                                     + std::pow(gp->supportNodes[i]->getY() -
                                                gp->supportNodes[j]->getY(), 2)
//                                      +std::pow(gp->supportNodes[i]->getZ() -
//                                                gp->supportNodes[j]->getZ(), 2)
                                     + std::pow((alpha_c * d_c), 2), q);
            }
        }
        break;

    case 1:
        // rbfType = 1 -> EXP
        // Compute R(i,j) = exp( -alpha_c*( (x_i - x_j)^2 + (y_i - y_j)^2 )
        //                / d_c^2 )
        qq = 1. / std::pow(alpha_c * d_c * 0.4, 2);
        for (auto i = 0u; i < nn; ++i) {
            for (auto j = 0u; j < nn; ++j) {
                g_o(i, j) = std::exp(-(std::pow(gp->supportNodes[i]->getX()
                                                - gp->supportNodes[j]->getX(), 2)
                                       + std::pow(gp->supportNodes[i]->getY()
                                                  - gp->supportNodes[j]->getY(), 2)
//                                        + std::pow (gp->supportNodes[i]->getZ()
//                                                  - gp->supportNodes[j]->getZ() , 2)
                                      ) * qq);
            }
        }
        break;

    case 10:
        // rbfType = 10 -> Wu-C2
        // Compute R(i,j) = (1 - radius/d_c)^5 * (8 + 40*(radius/d_c) + 48*(radius/d_c)^2 + 25*(radius/d_c)^3 + 5*(radius/d_c)^4 )
        for (auto i = 0u; i < nn; ++i) {
            for (auto j = 0u; j < nn; ++j) {
                radius = std::pow
                        (std::pow(gp->supportNodes[i]->getX() - gp->supportNodes[j]->getX(), 2)
                         + std::pow(gp->supportNodes[i]->getY() - gp->supportNodes[j]->getY(), 2)
//                            +std::pow( gp->supportNodes[i]->getZ() - gp->supportNodes[j]->getZ(), 2)
                        , 0.5);
                if (radius < d_c) {
                    g_o(i, j) = std::pow(1 - radius / d_c, 5)
                                * (8 +
                                   40 * (radius / d_c) +
                                   48 * std::pow(radius / d_c, 2) +
                                   25 * std::pow(radius / d_c, 3) +
                                   5 * std::pow(radius / d_c, 4));
                } else {
                    g_o(i, j) = 0.;
                }
            }
        }
        break;

    default:
        throw std::logic_error("Invalid value for rbfType");
    }
    if (mm > 0) {
        for (auto i = 0u; i < nn; ++i) {
            g_o(i, nn) = 1.;
            if (mm > 1) {
                g_o(i, nn + 1) = gp->supportNodes[i]->getX();
                g_o(nn + 1, i) = g_o(i, nn + 1);
                if (mm > 2) {
                    g_o(i, nn + 1) = gp->supportNodes[i]->getY();
                    g_o(nn + 1, i) = g_o(i, nn + 1);
                }
                else { std::cout << "WARNING: incomplete polynomial." << std::endl; }
            }
        }
    }
}

void ShapeFunctionRBF::computePhi(double xp, double yp, double zp)
{
    lmx::Vector<double> rhs(nn);
    lmx::Vector<double> phii(nn);
    lmx::DenseMatrix<double> rk;
    if (dim == 2) {
        rk.resize(6, nn); /**< Radial basis and derivatives. */
    } else if (dim == 3) {
        rk.resize(4, nn);
    } /**< Radial basis and derivatives. */
    double rr, rr2, xn, yn, zn;

    // For each derivative, first the rhs is calculated (= d^n/dx^n{ R(x) Pm(x) }^T) and then the linear system is solved to get {phi(n)}
    if (dim != 2)
        std::cout << "Warning, space dimension is partially implemented."
                  << std::endl;
    // Switch type of RBF:
    switch (rbfType) {
    case 0:
        // rbfType = 0 -> MQ
        // Compute R(i,j) = ((x_i - x_j)+(y_i - y_j) + (alpha_c*d_c)^2 )^q
        for (auto j = 0u; j < nn; ++j) {
            xn = gp->supportNodes[j]->getX();
            yn = gp->supportNodes[j]->getY();
            zn = 0;
//             zn = gp->supportNodes[j]->getZ();
            rr2 = std::pow(xp - xn, 2) + std::pow(yp - yn, 2) + std::pow(zp - zn, 2);
            rk(0, j) = std::pow(rr2 + std::pow(alpha_c * d_c, 2), q);
            rk(1, j) = 2. * q * std::pow(rr2 + std::pow(alpha_c * d_c, 2), q - 1.) * (xp - xn);
            rk(2, j) = 2. * q * std::pow(rr2 + std::pow(alpha_c * d_c, 2), q - 1.) * (yp - yn);
            if (dim == 2) {
                rk(3, j) = 2. * q * std::pow(rr2 + std::pow(alpha_c * d_c, 2), q - 1) +
                           4 * q * (q - 1) * std::pow(rr2 + std::pow(alpha_c * d_c, 2), q - 2) * std::pow(xp - xn, 2);
                // TODO: Check sign of derivatives
                rk(4, j) = -4 * q * (q - 1) * (xp - xn) * (yp - yn) *
                           std::pow(rr2 + std::pow(alpha_c * d_c, 2), q - 2);
                rk(5, j) = 2 * q * std::pow(rr2 + std::pow(alpha_c * d_c, 2), q - 1) +
                           4 * q * (q - 1) * std::pow(rr2 + std::pow(alpha_c * d_c, 2), q - 2) * std::pow(yp - yn, 2);
            }
            else if (dim == 3)
                rk(3, j) = 2. * q * std::pow(rr2 + std::pow(alpha_c * d_c, 2), q - 1.) * (zp - zn);
        }
        break;

    case 1: {
        // rbfType = 1 -> EXP
        // Compute R(i,j) = exp( -alpha_c*( (x_i - x_j)^2+(y_i - y_j)^2 ) / d_c^2 )
//      double qq = alpha_c / std::pow(d_c, 2);
        double qq = 1. / std::pow(alpha_c * d_c * 0.4, 2);
        for (auto j = 0u; j < nn; ++j) {
            xn = gp->supportNodes[j]->getX();
            yn = gp->supportNodes[j]->getY();
            zn = gp->supportNodes[j]->getZ();
            rr2 = std::pow(xp - xn, 2) + std::pow(yp - yn, 2) + std::pow(zp - zn, 2);
            rk(0, j) = std::exp(-qq * rr2);
            rk(1, j) = 2 * qq * std::exp(-qq * rr2) * (xp - xn);
            rk(2, j) = 2 * qq * std::exp(-qq * rr2) * (yp - yn);
            if (dim == 2) {
                rk(3, j) = 2 * qq * std::exp(-qq * rr2)
                           + 4 * qq * qq * std::pow(xp - xn, 2) * std::exp(-qq * rr2);
                rk(4, j) = -4 * qq * qq * std::exp(-qq * rr2) * (xp - xn) * (yp - yn);
                rk(5, j) = 2 * qq * std::exp(-qq * rr2)
                           + 4 * qq * qq * std::pow(yp - yn, 2) * std::exp(-qq * rr2);
            }
            else if (dim == 3) {
                rk(3, j) = 2 * qq * std::exp(-qq * rr2) * (zp - zn);
            }
        }
        break;
    }

    case 10:
        // rbfType = 10 -> Wu-C2
        // Compute R(i,j) = (1 - radius/d_c)^5 * (8 + 40*(radius/d_c) + 48*(radius/d_c)^2 + 25*(radius/d_c)^3 + 5*(radius/d_c)^4 )
        for (auto j = 0u; j < nn; ++j) {
            xn = gp->supportNodes[j]->getX();
            yn = gp->supportNodes[j]->getY();
            zn = gp->supportNodes[j]->getZ();
            rr2 = std::pow(xp - xn, 2) + std::pow(yp - yn, 2) + std::pow(zp - zn, 2);
            rr = std::pow(rr2, 0.5);

            if (rr < d_c) {
                rk(0, j) = std::pow(1 - rr / d_c, 5)
                           * (8 +
                              40 * (rr / d_c) +
                              48 * std::pow(rr / d_c, 2) +
                              25 * std::pow(rr / d_c, 3) +
                              5 * std::pow(rr / d_c, 4));
                rk(1, j) = (5 * std::pow(1 - rr / d_c, 4) * (-1 / d_c)
                            * (8 +
                               40 * (rr / d_c) +
                               48 * std::pow(rr / d_c, 2) +
                               25 * std::pow(rr / d_c, 3) +
                               5 * std::pow(rr / d_c, 4))
                            + std::pow(1 - rr / d_c, 5)
                              * (40 / d_c +
                                 96 * rr / std::pow(d_c, 2) +
                                 75 * std::pow(rr, 2) / std::pow(d_c, 3) +
                                 20 * std::pow(rr, 3) / std::pow(d_c, 4)))
                           * (xp - xn) / rr;
                rk(2, j) = (5 * std::pow(1 - rr / d_c, 4) * (-1 / d_c)
                            * (8 +
                               40 * (rr / d_c) +
                               48 * std::pow(rr / d_c, 2) +
                               25 * std::pow(rr / d_c, 3) +
                               5 * std::pow(rr / d_c, 4))
                            + std::pow(1 - rr / d_c, 5)
                              * (40 / d_c +
                                 96 * rr / std::pow(d_c, 2) +
                                 75 * std::pow(rr, 2) / std::pow(d_c, 3) +
                                 20 * std::pow(rr, 3) / std::pow(d_c, 4)))
                           * (yp - yn) / rr;
                if (dim == 3) {
                    rk(2, j) = (5 * std::pow(1 - rr / d_c, 4) * (-1 / d_c)
                                * (8 +
                                   40 * (rr / d_c) +
                                   48 * std::pow(rr / d_c, 2) +
                                   25 * std::pow(rr / d_c, 3) +
                                   5 * std::pow(rr / d_c, 4))
                                + std::pow(1 - rr / d_c, 5)
                                  * (40 / d_c +
                                     96 * rr / std::pow(d_c, 2) +
                                     75 * std::pow(rr, 2) / std::pow(d_c, 3) +
                                     20 * std::pow(rr, 3) / std::pow(d_c, 4)))
                               * (zp - zn) / rr;
                }
            } else {
                rk(0, j) = 0.;
                rk(1, j) = 0.;
                rk(2, j) = 0.;
                rk(3, j) = 0.;
                rk(4, j) = 0.;
                rk(5, j) = 0.;
            }
        }
        break;

    default:
        throw std::logic_error("Invalid value for rbfType");
    }
    if (mm > 0) {
        rk(0, nn) = 1.;
        if (mm > 1) {
            rk(0, nn + 1) = xp;
            rk(1, nn + 1) = 1.;
            if (mm > 2) {
                rk(0, nn + 2) = yp;
                rk(2, nn + 2) = 1.;
            }
            else std::cout << "WARNING: incomplete polynomial." << std::endl;
        }
    }

    // Solving linear system g_o*phii = rk(i) = rhs;
    lmx::LinearSystem<double> lsys(g_o, phii, rhs);
    // PHI and derivatives (6 systems are solved):
    //////////////////////////////////////
    // CHANGED FOR SOLVING ONLY dim+1 SYSTEMS: PHI AND FIRST PARTIAL
    // DERIVATIVES (dX, dY and dZ in 3D)
    //////////////////////////////////////
    auto n_sys = dim + 1;
    for (auto k = 0u; k < n_sys; ++k) {
        for (auto j = 0u; j < nn; ++j)
            rhs(j) = rk(k, j);
        lsys.solveYourself();
        for (auto j = 0u; j < nn; ++j) {
            this->phi(k, j) = phii(j);
        }
    }
    double sumphi = 0;
    for (auto i = 0u; i < nn; ++i) {
        sumphi += this->phi(0, i);
    }
    cout << "nPs: " << nn << ", SUM_PHI = " << sumphi << endl;
    cout << "PHI = " << phi << endl;


//  this->outputValues();
//  this->gnuplotOut();
}

} //Namespace mknix
