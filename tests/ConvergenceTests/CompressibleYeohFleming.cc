/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2010 - 2019 by the CMFE authors
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



#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/physics/elasticity/standard_tensors.h>
#include <deal.II/physics/elasticity/kinematics.h>

#include <iostream>
#include <string>

#include "ConvergenceTest.hpp"



using namespace dealii;


int main()
{
	const unsigned int dim = 3;

	const std::string material_parameter_file(std::string(SOURCE_DIR) + "/CompressibleYeohFleming.prm");

	Tensor<2,dim> tensor_F_0;

	tensor_F_0[0][0] = 0.93; tensor_F_0[1][1] =0.98 ; tensor_F_0[2][2] = 0.92;
/*		tensor_F_0[0][1] = 0; tensor_F_0[0][2] = 0;
		tensor_F_0[1][0] = 0; tensor_F_0[1][2] = 0;
		tensor_F_0[2][0] = 0; tensor_F_0[2][1] = 0;*/

	Tensor<2,dim> tensor_F_inf;

	tensor_F_inf[0][0] =0.95; tensor_F_inf[1][1] = 0.95; tensor_F_inf[2][2] = 0.95;
	tensor_F_inf[0][1] = 0.04; tensor_F_inf[0][2] = -0.07;
	tensor_F_inf[1][0] = -0.03; tensor_F_inf[1][2] = -0.09;
	tensor_F_inf[2][0] = 0; tensor_F_inf[2][1] = -0.1;

	const SymmetricTensor<2,dim> tensor_C_0 = Physics::Elasticity::Kinematics::C<dim>(tensor_F_0);

	const SymmetricTensor<2,dim> tensor_C_inf = Physics::Elasticity::Kinematics::C<dim>(tensor_F_inf);

	ConstitutiveLaws::convergence_test(material_parameter_file,
									   tensor_C_0,
									   tensor_C_inf);

	return 0;
}
