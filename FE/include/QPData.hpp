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

#pragma once


#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/lac/vector.h>

#include <boost/serialization/access.hpp>


#include <algorithm>


#include "BaseMaterial.hpp"


using namespace dealii;


namespace FE
{

class MaterialData
{
public:

	MaterialData(const std::string &material_parameter_file);

	MaterialData(const MaterialData &cpy);

	MaterialData& operator=(const MaterialData &cpy);

	std::string material_name;

	static void declare_parameters(ParameterHandler &prm);

	void parse_parameters(ParameterHandler &prm);

	// material parameters for Neo-Hookean material model
	double lambda;

	double mu;

	// material parameters for Ogden material model

	Vector<double> alpha;

	Vector<double> mu_O;

	int size;

	// material parameter Yeoh-Fleming model

	double A;

	double B;

	double C;

	double I_m;

	//material parameter incompressible Material

	double alpha_vol;

	double B_vol;

	template <class Archive>
	void serialize(Archive &ar,
	   	   	  	   const unsigned int version);
};


template <int dim>
class QPData
{
private:

	std::shared_ptr<ConstitutiveLaws::BaseMaterial<dim>> material;

public:

	QPData();

	QPData(const MaterialData &material_data);

	virtual ~QPData() = default;

	void reinit(const MaterialData &material_data);

	bool update_values(const Tensor<2,dim> &disp_grad);

	bool update_values(const SymmetricTensor<2,dim> &tensor_C);

	Tensor<2,dim> tensor_F;

	double det_J=0;

	SymmetricTensor<2,dim> tensor_C;

	SymmetricTensor<2,dim> tensor_S, tensor_S_iso, tensor_S_vol;

	SymmetricTensor<4,dim> tangent, tangent_iso, tangent_vol;
};

}//namespace FE
