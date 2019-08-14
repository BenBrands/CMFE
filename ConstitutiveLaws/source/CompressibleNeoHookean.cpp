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


#include <memory>


#include "CompressibleNeoHookean.hpp"


using namespace dealii;


namespace ConstitutiveLaws
{

template <int dim>
CompressibleNeoHookean<dim>::CompressibleNeoHookean(const double lambda_,
													const double mu_)
:
lambda(lambda_),
mu(mu_)
{}


template <int dim>
void
CompressibleNeoHookean<dim>::stress_S(SymmetricTensor<2,dim> &tensor_S,
									  const SymmetricTensor<2,dim> &tensor_C) const
{
	(void)tensor_S;
	(void)tensor_C;
/*
 *
 *
 *
 *
 *
 *
 */
}


template <int dim>
void
CompressibleNeoHookean<dim>::material_tangent(SymmetricTensor<4,dim> &tangent,
									  	  	  const SymmetricTensor<2,dim> &tensor_C) const
{
	(void)tangent;
	(void)tensor_C;
/*
 *
 *
 *
 *
 *
 *
 */
}


template <int dim>
void
CompressibleNeoHookean<dim>::stress_S_iso(SymmetricTensor<2,dim> &tensor_S,
		  	  	  	  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	(void)tensor_S;
	(void)tensor_C;

	Assert(false,ExcMessage("The function stress_S_iso is not callable for CompressibleNeoHookean"));
}


template <int dim>
void
CompressibleNeoHookean<dim>::stress_S_vol(SymmetricTensor<2,dim> &tensor_S,
		  	  	  	  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	(void)tensor_S;
	(void)tensor_C;

	Assert(false,ExcMessage("The function stress_S_vol is not callable for CompressibleNeoHookean"));
}


template <int dim>
void
CompressibleNeoHookean<dim>::material_tangent_iso(SymmetricTensor<4,dim> &tangent,
							  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	(void)tangent;
	(void)tensor_C;

	Assert(false,ExcMessage("The function material_tangent_iso is not callable for CompressibleNeoHookean"));
}


template <int dim>
void
CompressibleNeoHookean<dim>::material_tangent_vol(SymmetricTensor<4,dim> &tangent,
							  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	(void)tangent;
	(void)tensor_C;

	Assert(false,ExcMessage("The function material_tangent_vol is not callable for CompressibleNeoHookean"));
}


}//namespace ConstitutiveLaws


namespace boost
{
	namespace serialization
	{
		template<class Archive, int dim>
		inline void save_construct_data(Archive & ar,
										const ConstitutiveLaws::CompressibleNeoHookean<dim> *object_NH,
										const unsigned int /*version*/)
		{
			// save data required to construct instance
			ar << object_NH->lambda;
			ar << object_NH->mu;
		}


		template<class Archive, int dim>
		inline void load_construct_data(Archive & ar,
										ConstitutiveLaws::CompressibleNeoHookean<dim> *object_NH,
										const unsigned int /*version*/)
		{
			// retrieve data from archive required to construct new instance
			double lambda=0, mu=0;
			ar >> lambda;
			ar >> mu;

			// invoke inplace constructor to initialize instance of class
			::new(object_NH)ConstitutiveLaws::CompressibleNeoHookean<dim>(lambda,mu);
		}
	}//namespace serialization
}//namespace boost


#include "CompressibleNeoHookean.inst"
