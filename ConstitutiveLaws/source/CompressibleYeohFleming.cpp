




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
 * Authors: Laura Ruhland, Universität Erlangen-Nürnberg, 2019
 */


#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/exceptions.h>


#include <memory>


#include "CompressibleYeohFleming.hpp"


using namespace dealii;


namespace ConstitutiveLaws
{

template <int dim>
SymmetricTensor<4,dim> tensorproduct(SymmetricTensor<2,dim> &A,SymmetricTensor<2,dim> &B)
{
	SymmetricTensor<4,dim> C;

	for (unsigned int i=0; i<dim; i++)
		{
		for (unsigned int j=0; j<dim; j++)
			{
			for (unsigned int k=0; k<dim; k++)
				{
				for (unsigned int l=0; l<dim; l++)
					{
						C[i][j][k][l]=(A[i][k]*B[j][l]+A[i][l]*B[j][k]);
					}
				}
			}
		}

	C*=  0.5;

	return C;
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

	SymmetricTensor<4,dim> tensor_P_bar = tensorproduct(inv_C,inv_C);

	tensor_P_bar -= (double)1/dim * outer_product(inv_C,inv_C);

	return tensor_P_bar;
}



template <int dim>
CompressibleYeohFleming<dim>::CompressibleYeohFleming(const double A_,
													const double B_,const double C_, const double I_m_, const double alpha_, const double B_vol_)
:
A(A_),
B(B_),
C(C_),
I_m(I_m_),
alpha(alpha_),
B_vol(B_vol_)
{}



template<int dim>
SymmetricTensor<2,dim>
CompressibleYeohFleming<dim>::stress_S_bar(const SymmetricTensor<2,dim> &tensor_C) const
{
	double J= sqrt(determinant(tensor_C));

	double I1=first_invariant(tensor_C);

	double I1_ = std::pow(J,-(double)2/dim) * I1;

	double factor=2*A*std::exp(-B*(I1_-(double)dim));
	factor += 2*C*(I_m-dim)/(I_m-I1_);

	SymmetricTensor<2,dim> tensor_S_bar = unit_symmetric_tensor<dim>();

	tensor_S_bar*=factor;

	return tensor_S_bar;
}



template<int dim>
SymmetricTensor<4,dim>
CompressibleYeohFleming<dim>::tangente_C_bar(const SymmetricTensor<2,dim> &tensor_C)const
{
	double J= sqrt(determinant(tensor_C));

	double I1=first_invariant(tensor_C);

	double I1_ = std::pow(J,-(double)2/dim) * I1;

	double factor = 4*std::pow(J,-(double)4/dim)*(-A*B*std::exp(-B*(I1_-dim)));

	factor += 4*std::pow(J,-(double)4/dim)*(C*(I_m-dim)/std::pow(I_m-I1_,2));

	SymmetricTensor<4,dim> tensor_C_bar = outer_product(unit_symmetric_tensor<dim>(),
			  unit_symmetric_tensor<dim>());

	tensor_C_bar*=factor;

	return tensor_C_bar;
}



template <int dim>
void
CompressibleYeohFleming<dim>::stress_S(SymmetricTensor<2,dim> &tensor_S,
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
CompressibleYeohFleming<dim>::material_tangent(SymmetricTensor<4,dim> &tangent,
									  	  	  const SymmetricTensor<2,dim> &tensor_C) const
{
	SymmetricTensor<4,dim> tangent_iso;

	SymmetricTensor<4,dim> tangent_vol;

	material_tangent_iso(tangent_iso,tensor_C);

	material_tangent_vol(tangent_vol,tensor_C);

	tangent = tangent_iso;

	tangent += tangent_vol;
}



template <int dim>
void
CompressibleYeohFleming<dim>::stress_S_iso(SymmetricTensor<2,dim> &tensor_S,
		  	  	  	  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	double J=sqrt(determinant(tensor_C));

	SymmetricTensor<2,dim> tensor_S_bar = stress_S_bar(tensor_C);

	tensor_S=projection_tensor_P(tensor_C)*tensor_S_bar;

	tensor_S*= std::pow(J,-(double)2/dim);
}



template <int dim>
void
CompressibleYeohFleming<dim>::stress_S_vol(SymmetricTensor<2,dim> &tensor_S,
		  	  	  	  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	double J= sqrt(determinant(tensor_C));

	tensor_S = invert(tensor_C);

	tensor_S *= J*(B_vol/alpha)*std::sinh(alpha*(J-1));
}



template <int dim>
void
CompressibleYeohFleming<dim>::material_tangent_iso(SymmetricTensor<4,dim> &tangent,
							  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	double J;

	J=sqrt(determinant(tensor_C));

	SymmetricTensor<2,dim> inv_C= invert(tensor_C);

	SymmetricTensor<2,dim> S_iso;

	stress_S_iso(S_iso,tensor_C);

	SymmetricTensor<2,dim> S_bar = stress_S_bar(tensor_C);

	SymmetricTensor<4,dim> C_bar = tangente_C_bar(tensor_C);

	SymmetricTensor<4,dim> proj_tensor_P = projection_tensor_P(tensor_C);

	SymmetricTensor<4,dim> mod_proj_tensor_P = modified_projection_tensor_P(tensor_C);

	tangent=proj_tensor_P*C_bar*transpose(proj_tensor_P);

	const double factor = (double)2/dim * std::pow(J,-(double)2/dim)*(S_bar*tensor_C);

	tangent += mod_proj_tensor_P*factor;

	tangent -= (double)2/dim * (outer_product(inv_C,S_iso)+outer_product(S_iso,inv_C));
}


template <int dim>
void
CompressibleYeohFleming<dim>::material_tangent_vol(SymmetricTensor<4,dim> &tangent,
							  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	double J = sqrt(determinant(tensor_C));

	const double sinh_a_J = std::sinh(alpha*(J-1));

	const double cosh_a_J = std::cosh(alpha*(J-1));

	SymmetricTensor<4,dim> C_vol1;
	SymmetricTensor<4,dim> C_vol2;

	SymmetricTensor<2,dim> C_inv= invert(tensor_C);

	C_vol1=outer_product(C_inv,C_inv);

	double factor = (B_vol/alpha)*J*sinh_a_J;

	factor += std::pow(J,2)*B_vol*cosh_a_J;

	C_vol1 *= factor;

	C_vol2=tensorproduct(C_inv,C_inv);

	C_vol2*=2*J*(B_vol/alpha)*sinh_a_J;

	tangent=C_vol1-C_vol2;
}



}//namespace ConstitutiveLaws


namespace boost
{
	namespace serialization
	{
		template<class Archive, int dim>
		inline void save_construct_data(Archive & ar,
										const ConstitutiveLaws::CompressibleYeohFleming<dim> *object_YF,
										const unsigned int /*version*/)
		{
			// save data required to construct instance
			ar << object_YF->A;
			ar << object_YF->B;
			ar << object_YF->C;
			ar << object_YF->I_m;
			ar << object_YF->alpha;
			ar << object_YF->B_vol;
		}


		template<class Archive, int dim>
		inline void load_construct_data(Archive & ar,
										ConstitutiveLaws::CompressibleYeohFleming<dim> *object_YF,
										const unsigned int /*version*/)
		{
			// retrieve data from archive required to construct new instance
			double A=0, B=0, C=0, I_m=0, alpha=0, B_vol=0;
			ar >> A;
			ar >> B;
			ar >> C;
			ar >> I_m;
			ar >> alpha;
			ar >> B_vol;

			// invoke inplace constructor to initialize instance of class
			::new(object_YF)ConstitutiveLaws::CompressibleYeohFleming<dim>(A,B,C,I_m,alpha,B_vol);
		}
	}//namespace serialization
}//namespace boost


#include "CompressibleYeohFleming.inst"
