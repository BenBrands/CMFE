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
#include <deal.II/physics/elasticity/standard_tensors.h>


#include <memory>


#include "CompressibleNeoHookean.hpp"


using namespace dealii;



template <int dim>
SymmetricTensor<4,dim> dInvC_dC(const SymmetricTensor<2,dim> &inv_C)
{
	SymmetricTensor<4,dim> result;

	for (unsigned int i=0; i<dim; ++i)
		for (unsigned int j=i; j<dim; ++j)
			for (unsigned int k=0; k<dim; ++k)
				for (unsigned int l=k; l<dim; ++l)
					result[i][j][k][l] = (inv_C[i][k] * inv_C[j][l] + inv_C[i][l] * inv_C[j][k]) / 2;

	return result;
}



template <int dim>
SymmetricTensor<4,dim> projection_tensor_P(const SymmetricTensor<2,dim> &tensor_C)
{
	SymmetricTensor<2,dim> inv_C = invert(tensor_C);

	SymmetricTensor<4,dim> tensor_P = outer_product(inv_C,tensor_C);

	tensor_P *= -(double)1/dim;

	for (unsigned int i=0; i<dim; ++i)
		for (unsigned int j=0; j<dim; ++j)
			for (unsigned int k=0; k<dim; ++k)
				for (unsigned int l=0; l<dim; ++l)
					if (i==k && j==l)
						tensor_P[i][j][k][l] += 1;

	return tensor_P;
}



template <int dim>
SymmetricTensor<4,dim> modified_projection_tensor_P(const SymmetricTensor<2,dim> &tensor_C)
{
	SymmetricTensor<2,dim> inv_C = invert(tensor_C);

	SymmetricTensor<4,dim> tensor_P_bar = dInvC_dC(inv_C);

	tensor_P_bar -= (double)1/dim * outer_product(inv_C,inv_C);

	return tensor_P_bar;
}



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
	SymmetricTensor<2,dim> S_iso;

	SymmetricTensor<2,dim> S_vol;

	stress_S_iso(S_iso,tensor_C);

	stress_S_vol(S_vol,tensor_C);

	tensor_S = S_iso;

	tensor_S += S_vol;
}



template <int dim>
void
CompressibleNeoHookean<dim>::material_tangent(SymmetricTensor<4,dim> &tangent,
									  	  	  const SymmetricTensor<2,dim> &tensor_C) const
{
	tangent=0;

	SymmetricTensor<4,dim> tangent_iso;

	SymmetricTensor<4,dim> tangent_vol;

	material_tangent_iso(tangent_iso,tensor_C);

	material_tangent_vol(tangent_vol,tensor_C);

	tangent = tangent_iso;

	tangent += tangent_vol;
}



template <int dim>
void
CompressibleNeoHookean<dim>::stress_S_iso(SymmetricTensor<2,dim> &tensor_S,
		  	  	  	  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	double J;

	J=sqrt(determinant(tensor_C));

	SymmetricTensor<2,dim> tensor_S_bar = mu*unit_symmetric_tensor<dim>();

	tensor_S=projection_tensor_P(tensor_C)*tensor_S_bar;

	tensor_S*= std::pow(J,-(double)2/dim);
}



template <int dim>
void
CompressibleNeoHookean<dim>::stress_S_vol(SymmetricTensor<2,dim> &tensor_S,
		  	  	  	  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	 double J= sqrt(determinant(tensor_C));

	tensor_S = invert(tensor_C);

	double kappa=lambda+(double)2/3*mu;

	double dpsi_vol_dJ=kappa/2*(J-std::pow(J,-1.0));

	tensor_S *= J*dpsi_vol_dJ;
}



template <int dim>
void
CompressibleNeoHookean<dim>::material_tangent_iso(SymmetricTensor<4,dim> &tangent,
							  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	tangent=0;

	double J;

	J=sqrt(determinant(tensor_C));

	SymmetricTensor<2,dim> inv_C= invert(tensor_C);

	SymmetricTensor<2,dim> S_bar =mu*unit_symmetric_tensor<dim>();

	SymmetricTensor<2,dim> S_iso=projection_tensor_P(tensor_C)*S_bar;

	S_iso*= std::pow(J,-(double)2/dim);

	SymmetricTensor<2,dim> tensor_zero;

	SymmetricTensor<4,dim> C_bar = outer_product(tensor_zero,tensor_zero);

	SymmetricTensor<4,dim> proj_tensor_P = projection_tensor_P(tensor_C);

	SymmetricTensor<4,dim> mod_proj_tensor_P = modified_projection_tensor_P(tensor_C);

	tangent=proj_tensor_P*C_bar*transpose(proj_tensor_P);

	double factor = ((double)2/dim) * std::pow(J,-(double)2/dim)*(S_bar*tensor_C);

	tangent += mod_proj_tensor_P*factor;

	tangent -= ((double)2/dim) * (outer_product(inv_C,S_iso)+outer_product(S_iso,inv_C));
}



template <int dim>
void
CompressibleNeoHookean<dim>::material_tangent_vol(SymmetricTensor<4,dim> &tangent,
							  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	 double J = sqrt(determinant(tensor_C));

	SymmetricTensor<4,dim> C_vol1;
	SymmetricTensor<4,dim> C_vol2;

	SymmetricTensor<2,dim> C_inv= invert(tensor_C);

	C_vol1=outer_product(C_inv,C_inv);

	double kappa=lambda+((double)2/3)*mu;

	double factor =kappa*J*J;

	C_vol1 *= factor;

	C_vol2=dInvC_dC(C_inv);

	C_vol2*=kappa*(1-J*J);

	tangent=C_vol1+C_vol2;
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
