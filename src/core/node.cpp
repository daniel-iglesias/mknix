/******************************************************************************
 *  Copyright (C) 2015 by Daniel Iglesias                                     *
 *                                                                            *
 *  This file is part of MkniX.                                             *
 *                                                                            *
 *  MkniX is free software: you can redistribute it and/or modify           *
 *  it under the terms of the GNU Lesser General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  MkniX is distributed in the hope that it will be useful,                *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Lesser General Public License for more details.                       *
 *                                                                            *
 *  You should have received a copy of the GNU Lesser General Public          *
 *  License along with MkniX.  If not, see <http://www.gnu.org/licenses/>.  *
 *****************************************************************************/

#include "node.h"
#include "shapefunction.h"

namespace mknix
{

/**
 * @brief Default constructor for Node.
 */
Node::Node() { }

/**
 * @brief Copy constructor. Copies position, displacement, and temperature from another Node.
 * @param p_node_in Reference to the source Node.
 */
Node::Node(const Node& p_node_in)
    : Point(p_node_in)
    , qx(p_node_in.getqx(0))
    , qy(p_node_in.getqx(1))
    , qz(p_node_in.getqx(2))
    , qt(p_node_in.getqt())
    , weight(0) { }

/**
 * @brief Pointer-based copy constructor. Copies position, displacement, and temperature from another Node.
 * @param p_node_in Pointer to the source Node.
 */
Node::Node(const Node* p_node_in)
    : Point(p_node_in)
    , qx(p_node_in->getqx(0))
    , qy(p_node_in->getqx(1))
    , qz(p_node_in->getqx(2))
    , qt(p_node_in->getqt())
    , weight(0) { }

/**
 * @brief Constructs a Node at a given position with initial displacement equal to its coordinates.
 * @param i_in Global node index.
 * @param x_in X coordinate.
 * @param y_in Y coordinate.
 * @param z_in Z coordinate.
 */
Node::Node(int i_in, double x_in, double y_in, double z_in)
    : Point(i_in, x_in, y_in, z_in)
    , qx(x_in)
    , qy(y_in)
    , qz(z_in)
    , qt(0)
    , weight(0)
{
}

/**
 * @brief Destructor for Node.
 */
Node::~Node()
{
    /*std::cout << "---DESTROYED NODE---" << std::endl;*/
}

/**
 * @brief Returns the current configuration coordinate for the given degree of freedom.
 *        For meshfree nodes (MLS, 1D, 2D, 3D), evaluates the interpolated value from support nodes;
 *        for FEM nodes, returns the stored displacement directly.
 * @param gdl Degree of freedom index (0 = x, 1 = y, 2 = z).
 * @return Interpolated or direct configuration value.
 */
double Node::getConf(int gdl) const
{
    // if not delta_kronecker : x = sum_i( phi_i * q_i )
    double conf_value(0);
//   cout << this->num << endl;
    if (shapeFunType == "MLS"
            || shapeFunType == "1D"
            || shapeFunType == "2D"
            || shapeFunType == "3D")
    {
        for (auto i = 0u; i < supportNodesSize; ++i)
        {
            conf_value += shapeFun->getPhi(0, i) * (supportNodes[i]->getqx(gdl));
//       cout << endl << "NODE " << num << " conf += (" << shapeFun->getPhi(0, i) <<" * "
// 	   << (supportNodes[i]->getx(gdl) ) <<" ) = " << conf_value << endl;

        }
//    cout << "conf_value = " << conf_value << ", q = " << getx(gdl) << endl;
    }
    else
    {
        conf_value = this->getqx(gdl);
    }
    return conf_value;
}

/**
 * @brief Returns the current temperature at this node.
 *        For MLS nodes, interpolates from support node temperatures;
 *        for FEM nodes, returns the stored temperature directly.
 * @return Temperature value.
 */
double Node::getTemp() const
{
    // if not delta_kronecker : x = sum_i( phi_i * q_i )
    double conf_value(0);
//   cout << this->num << endl;
    if (shapeFunType == "MLS")
    {
        for (auto i = 0u; i < supportNodesSize; ++i)
        {
            conf_value += shapeFun->getPhi(0, i) * (supportNodes[i]->getqt());
//       cout << endl << "conf += (" << shapeFun->getPhi(0, i) <<" * "
//         << (supportNodes[i]->getx(gdl) ) <<" ) = " << conf_value << endl;

        }
//    cout << "conf_value = " << conf_value << ", qt = " << getqt() << endl;
    }
    else
    {
        conf_value = this->getqt();
    }
    return conf_value;
}

/**
 * @brief Returns the number of support nodes for a given derivative order.
 * @param deriv Derivative order (0 = value, >0 = derivative).
 * @return Number of support nodes, or 1 for standard FEM nodes.
 */
size_t Node::getSupportSize(int deriv)
{
//     cout << this->shapeFunType << endl;
    if (deriv == 0)
    {
        // if not delta_kronecker : x = sum_i( phi_i * q_i )
        if (shapeFunType == "MLS"
                || shapeFunType == "1D"
                || shapeFunType == "2D"
                || shapeFunType == "3D"
           )
        {
            return this->supportNodes.size();
        }
        else
        {
            return 1;
        }
    }
    else   // derivative order > 1
    {
        if (shapeFunType == "MLS" || shapeFunType == "RBF" || shapeFunType == "1D")
        {
            return this->supportNodes.size();
        }
        else
        {
            return static_cast<size_t>(1E10); // Produce an infinite loop
        }
    }

}

/**
 * @brief Returns the global node number of a support node at a given derivative order and index.
 * @param deriv Derivative order.
 * @param s_node Index of the support node within the local support set.
 * @return Global node number of the support node.
 */
int Node::getSupportNodeNumber(int deriv, int s_node)
{
    if (deriv == 0)
    {
        // if not delta_kronecker : x = sum_i( phi_i * q_i )
        if (shapeFunType == "MLS"
                || shapeFunType == "1D"
                || shapeFunType == "2D"
                || shapeFunType == "3D")
        {
            return this->supportNodes[s_node]->getNumber();
        }
        else
        {
            return this->getNumber();
        }
    }
    else   // derivative order > 1
    {
        if (shapeFunType == "MLS"
                || shapeFunType == "RBF"
                || shapeFunType == "1D"
                || shapeFunType == "2D"
                || shapeFunType == "3D")
        {
            return this->supportNodes[s_node]->getNumber();
        }
        else
        {
            return -1; // Produce an error
        }
    }

}

/**
 * @brief Returns the shape function value (phi) for a given derivative order and support node index.
 * @param deriv Derivative order (0 = value, 1 = dx, 2 = dy, etc.).
 * @param s_node Index of the support node within the local support set.
 * @return Shape function value; 1.0 for standard FEM nodes at derivative order 0.
 */
double Node::getShapeFunValue(int deriv, int s_node)
{
    if (deriv == 0
            && (shapeFunType != "MLS"
                && shapeFunType != "1D"
                && shapeFunType != "2D"
                && shapeFunType != "3D"))
    {
        return 1.;
    }
    else
    {
        return this->shapeFun->getPhi(deriv, s_node);
    }
}


/**
 * @brief Updates the node's displacement components from a global configuration vector.
 * @param globalConf Global displacement state vector.
 * @param dim Spatial dimension (2 or 3).
 */
void Node::setqx(const lmx::Vector<data_type>& globalConf, int dim)
{
    qx = globalConf.readElement(dim * num);
    qy = globalConf.readElement(dim * num + 1);
    if (dim == 3)
    {
        qz = globalConf.readElement(dim * num + 2);
    }
}

/**
 * @brief Updates the node's temperature from a global temperature vector.
 * @param globalTemp Global temperature state vector.
 */
void Node::setqt(const lmx::Vector<data_type>& globalTemp)
{
    qt = globalTemp.readElement(thermalNum);
}

/* Function needed to initialize self-supported nodes that are part of
 * the formulation, but need to be defined at the reference configuration
 * in a different space from the Lagrangian (ie. RB director vectors) */
//   void Node::setShapeFunValue( double the_value, int deriv, int s_node )
//   {
//     this->shapeFun->setPhi( the_value, deriv, s_node );
//   }

}

