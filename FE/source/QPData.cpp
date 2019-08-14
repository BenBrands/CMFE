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


#include <deal.II/base/patterns.h>

#include <deal.II/physics/elasticity/standard_tensors.h>
#include <deal.II/physics/elasticity/kinematics.h>


#include "QPData.hpp"

#include "CompressibleNeoHookean.hpp"


using namespace dealii;


namespace FE
{

MaterialData::MaterialData(const std::string &material_parameter_file)
{
	ParameterHandler prm;

	declare_parameters(prm);

	prm.parse_input(material_parameter_file);

	parse_parameters(prm);
}


MaterialData::MaterialData(const MaterialData &cpy)
:
material_name(cpy.material_name),
lambda(cpy.lambda),
mu(cpy.mu)
{}


MaterialData&
MaterialData::operator=(const MaterialData &cpy)
{
	if (this != &cpy)
	{
		material_name = cpy.material_name;

		lambda = cpy.lambda;

		mu = cpy.mu;
	}
	return *this;
}


void
MaterialData::declare_parameters(ParameterHandler &prm)
{
	prm.enter_subsection("Type of Material");
	{
		Patterns::Selection selection("NeoHookean|Ogden|Yeoh|TNM|Gent|ArrudaBoyce");

	    prm.declare_entry("Name",
	                      "NeoHookean",
	                      selection,
	                      "identifier for material model");
	}
	prm.leave_subsection();

	prm.enter_subsection("Material properties NeoHookean");
	{
		prm.declare_entry("First Lame Parameter",
						  "12",
						  Patterns::Double(0.));

		prm.declare_entry("Second Lame Parameter",
						  "8",
						  Patterns::Double(0.));
	}
	prm.leave_subsection();
}


void
MaterialData::parse_parameters(ParameterHandler &prm)
{
	prm.enter_subsection("Type of Material");
	{
		material_name = prm.get("Name");
	}
	prm.leave_subsection();

	if (material_name.compare("NeoHookean") == 0)
	{
		prm.enter_subsection("Material properties NeoHookean");
		{
			lambda = prm.get_double("First Lame Parameter");

	        mu = prm.get_double("Second Lame Parameter");
		}
		prm.leave_subsection();
	}
	else if (material_name.compare("Ogden") == 0)
	{
		Assert(false,ExcMessage("Implementation missing"));
	}
	else if (material_name.compare("Yeoh") == 0)
	{
		Assert(false,ExcMessage("Implementation missing"));
	}
	else if (material_name.compare("TNM") == 0)
	{
		Assert(false,ExcMessage("Implementation missing"));
	}
	else if (material_name.compare("GENT") == 0)
	{
		Assert(false,ExcMessage("Implementation missing"));
	}
	else if (material_name.compare("ArrudaBoyce") == 0)
	{
		Assert(false,ExcMessage("Implementation missing"));
	}
}


template <class Archive>
void
MaterialData::serialize(Archive &ar,
   	   	  	   	   	   	const unsigned int /*version*/)
{
	ar & material_name;

	ar & lambda;

	ar & mu;
}


template <int dim>
QPData<dim>::QPData()
:
tensor_F(Physics::Elasticity::StandardTensors<dim>::I),
det_J(1),
tensor_C(Physics::Elasticity::StandardTensors<dim>::I)
{}


template <int dim>
QPData<dim>::QPData(const MaterialData &material_data)
:
tensor_F(Physics::Elasticity::StandardTensors<dim>::I),
det_J(1),
tensor_C(Physics::Elasticity::StandardTensors<dim>::I)
{
	if (material_data.material_name.compare("NeoHookean")==0)
	{
		material = std::make_shared<ConstitutiveLaws::CompressibleNeoHookean<dim>>(material_data.lambda,material_data.mu);
	}
	else
	{
		AssertThrow(false,ExcMessage("Other material models have to be added"));
	}
}


template <int dim>
void
QPData<dim>::reinit(const MaterialData &material_data)
{
	tensor_C = Physics::Elasticity::StandardTensors<dim>::I;

	det_J = 1;

	if (material_data.material_name.compare("NeoHookean")==0)
	{
		material = std::make_shared<ConstitutiveLaws::CompressibleNeoHookean<dim>>(material_data.lambda,material_data.mu);
	}
	else
	{
		AssertThrow(false,ExcMessage("Other material models have to be added"));
	}
}


template <int dim>
void
QPData<dim>::update_values(const Tensor<2,dim> &disp_grad)
{
	tensor_F = Physics::Elasticity::Kinematics::F(disp_grad);

	tensor_C = Physics::Elasticity::Kinematics::C(tensor_F);

	det_J = sqrt(determinant(tensor_C));

	AssertThrow(det_J > 0,ExcMessage("det_J <= 0"));

	material->stress_S(tensor_S,tensor_C);

	material->material_tangent(tangent,tensor_C);
}


}//namespace FE

#include "QPData.inst"
