/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2010 - 2019 by the CMFE
 *
 * This file is part of the CMFE library.
 *
 * The CMFE library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of CMFE.
 *
 * ---------------------------------------------------------------------

 *
 * Authors: Benjamin Brands, Universität Erlangen-Nürnberg, 2019
 */

#pragma once



#include <deal.II/base/utilities.h>
#include <deal.II/base/symmetric_tensor.h>

#include "BaseMaterial.hpp"
#include "QPData.hpp"

#include <iostream>



using namespace dealii;

namespace ConstitutiveLaws
{

template <int dim>
void convergence_test(const std::string &material_parameter_file,
					  const SymmetricTensor<2,dim> &C_0,
					  const SymmetricTensor<2,dim> &C_inf)
{
	FE::MaterialData material(material_parameter_file);

	FE::QPData<dim> qpd_inf(material);

	qpd_inf.update_values(C_inf);

	FE::QPData<dim> qpd(material);

	qpd.update_values(C_0);

	SymmetricTensor<2,dim> residuum = qpd_inf.tensor_S - qpd.tensor_S;

	const double norm_initial_residuum = residuum.norm();

	std::cout << "Norm of initial residuum: " << std::setprecision(5) << norm_initial_residuum << std::endl;

	for (unsigned int it=0; it<10; ++it)
	{
		SymmetricTensor<4,dim> inv = invert(0.5 * qpd.tangent);

		SymmetricTensor<2,dim> delta_C = inv * residuum;

		qpd.update_values(qpd.tensor_C+delta_C);

		residuum = qpd_inf.tensor_S - qpd.tensor_S;

		const double norm_residuum = residuum.norm();

		std::cout << "   Norm of residuum after iteration " << it+1 << ": "
				  << std::setprecision(5) << norm_residuum << std::endl;

		if (norm_residuum/norm_initial_residuum < 1e-6)
		{
			std::cout << "Obtained convergence after " << it+1 << " iterations" << std::endl;

			break;
		}
	}
}

}//namespace ConstitutiveLaws
