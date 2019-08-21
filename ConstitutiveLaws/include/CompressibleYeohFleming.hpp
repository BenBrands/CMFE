
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
 * Authors: Laura Ruhland, Universität Erlangen-Nürnberg, 2019
 */

#pragma once


#include <deal.II/base/symmetric_tensor.h>


#include <boost/serialization/access.hpp>


#include "BaseMaterial.hpp"


using namespace dealii;


namespace ConstitutiveLaws
{

template <int dim>
class CompressibleYeohFleming : public BaseMaterial<dim>
{
private:

	friend class boost::serialization::access;

	template <class Archive>
	void load(Archive &/*ar*/,
			  const unsigned int /*version*/) {}

	template <class Archive>
	void save(Archive &/*ar*/,
			  const unsigned int /*version*/) const {}

	BOOST_SERIALIZATION_SPLIT_MEMBER()

public:

	CompressibleYeohFleming(const double A,
						   const double B, const double C, const double I_m);

	virtual ~CompressibleYeohFleming() {}

	virtual void stress_S(SymmetricTensor<2,dim> &tensor_S,
						  const SymmetricTensor<2,dim> &tensor_C) const;

	virtual void stress_S_iso(SymmetricTensor<2,dim> &tensor_S,
			  	  	  	  	  const SymmetricTensor<2,dim> &tensor_C) const;

	virtual void stress_S_vol(SymmetricTensor<2,dim> &tensor_S,
			  	  	  	  	  const SymmetricTensor<2,dim> &tensor_C) const;

	virtual void material_tangent(SymmetricTensor<4,dim> &tangent,
								  const SymmetricTensor<2,dim> &tensor_C) const;

	virtual void material_tangent_iso(SymmetricTensor<4,dim> &tangent,
								  	  const SymmetricTensor<2,dim> &tensor_C) const;

	virtual void material_tangent_vol(SymmetricTensor<4,dim> &tangent,
								  	  const SymmetricTensor<2,dim> &tensor_C) const;

	const double A, B, C, I_m;
};

}//namespace ConstitutiveLaws
