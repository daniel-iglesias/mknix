/***************************************************************************
 *   Copyright (C) 2013 by Daniel Iglesias                                 *
 *   http://code.google.com/p/mknix                                        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "node.h"
#include "shapefunction.h"

namespace mknix
{

Node::Node() { }

Node::Node(const Node& p_node_in)
    : Point(p_node_in)
    , qx(p_node_in.getqx(0))
    , qy(p_node_in.getqx(1))
    , qz(p_node_in.getqx(2))
    , qt(p_node_in.getqt())
    , weight(0) { }

Node::Node(const Node* p_node_in)
    : Point(p_node_in)
    , qx(p_node_in->getqx(0))
    , qy(p_node_in->getqx(1))
    , qz(p_node_in->getqx(2))
    , qt(p_node_in->getqt())
    , weight(0) { }

Node::Node(int i_in, double x_in, double y_in, double z_in)
    : Point(i_in, x_in, y_in, z_in)
    , qx(x_in)
    , qy(y_in)
    , qz(z_in)
    , qt(0)
    , weight(0)
{
}

Node::~Node()
{
    /*std::cout << "---DESTROYED NODE---" << std::endl;*/
}

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


void Node::setqx(const lmx::Vector<data_type>& globalConf, int dim)
{
    qx = globalConf.readElement(dim * num);
    qy = globalConf.readElement(dim * num + 1);
    if (dim == 3)
    {
        qz = globalConf.readElement(dim * num + 2);
    }
}

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

