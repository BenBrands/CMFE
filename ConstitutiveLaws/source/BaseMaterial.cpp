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


#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/exceptions.h>


#include <boost/mpi.hpp>


#include <memory>


#include "BaseMaterial.hpp"


using namespace dealii;


namespace ConstitutiveLaws
{

template <int dim>
void
BaseMaterial<dim>::stress_S_iso(SymmetricTensor<2,dim> &tensor_S,
		  	  	  	  	  	  	const SymmetricTensor<2,dim> &tensor_C) const
{
	(void)tensor_S;
	(void)tensor_C;

	Assert(false,ExcMessage("The function stress_S_iso is not implemented for the derived class"));
}


template <int dim>
void
BaseMaterial<dim>::stress_S_vol(SymmetricTensor<2,dim> &tensor_S,
		  	  	  	  	  	  	const SymmetricTensor<2,dim> &tensor_C) const
{
	(void)tensor_S;
	(void)tensor_C;

	Assert(false,ExcMessage("The function stress_S_vol is not implemented for the derived class"));
}


template <int dim>
void
BaseMaterial<dim>::material_tangent_iso(SymmetricTensor<4,dim> &tangent,
							  	  	  	const SymmetricTensor<2,dim> &tensor_C) const
{
	(void)tangent;
	(void)tensor_C;

	Assert(false,ExcMessage("The function material_tangent_iso is not implemented for the derived class"));
}


template <int dim>
void
BaseMaterial<dim>::material_tangent_vol(SymmetricTensor<4,dim> &tangent,
							  	  	  	const SymmetricTensor<2,dim> &tensor_C) const
{
	(void)tangent;
	(void)tensor_C;

	Assert(false,ExcMessage("The function material_tangent_vol is not implemented for the derived class"));
}


template <int dim>
template <class Archive>
void
BaseMaterial<dim>::serialize(Archive &ar,
			   	   	   	   	 const unsigned int version)
{
	(void)ar;
	(void)version;
}

}//namespace ConstitutiveLaws

#include "BaseMaterial.inst"
