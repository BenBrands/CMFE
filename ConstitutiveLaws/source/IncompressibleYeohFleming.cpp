
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


#include "IncompressibleYeohFleming.hpp"


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
										C[i][j][k][l]=0;

										C[i][j][k][l]=(A[i][k]*B[j][l]+A[i][l]*B[j][k]);
									}
								}
							}
			}

			C*=  0.5;

			return C;
		}

template <int dim>
SymmetricTensor<4,dim> product(SymmetricTensor<2,dim> &A,SymmetricTensor<2,dim> &B)
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
										C[i][j][k][l]=A[i][j]*B[k][l];
									}
								}
							}
			}

			return C;
		}

template<typename T>
int kroneckerdelta(T &x,T &y)
{
	int kd;

	if (x==y)
	{
		kd= 1;
	}
	else
	{
		kd= 0;
	}

	return kd;
}

template <int dim>
IncompressibleYeohFleming<dim>::IncompressibleYeohFleming(const double A_,
													const double B_,const double C_, const double I_m_, const double alpha_, const double B_vol_)
:
A(A_),
B(B_),
C(C_),
I_m(I_m_),
alpha(alpha_),
B_vol(B_vol_)
{}


template <int dim>
void
IncompressibleYeohFleming<dim>::stress_S(SymmetricTensor<2,dim> &tensor_S,
									  const SymmetricTensor<2,dim> &tensor_C) const
{
	(void)tensor_S;
	(void)tensor_C;

	Assert(false,ExcMessage("The function stress_S_iso is not callable for IncompressibleYeohFleming"));


}


template <int dim>
void
IncompressibleYeohFleming<dim>::material_tangent(SymmetricTensor<4,dim> &tangent,
									  	  	  const SymmetricTensor<2,dim> &tensor_C) const
{
	(void)tangent;
	(void)tensor_C;

	Assert(false,ExcMessage("The function material_tangent_iso is not callable for IncompressibleYeohFleming"));

}


template <int dim>
void
IncompressibleYeohFleming<dim>::stress_S_iso(SymmetricTensor<2,dim> &tensor_S,
		  	  	  	  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{

	 double J=1;

	 J=sqrt(determinant(tensor_C));

	Assert(J>0,ExcMessage("Determinant J is smaller or equal to 0"));

	SymmetricTensor<2,dim> inv_C= invert(tensor_C);

	double I1=first_invariant(tensor_C);

	double I1_=std::pow(J,(-2/(double)dim))*I1;

	SymmetricTensor<2,dim> I;

	I=unit_symmetric_tensor<dim>();

	tensor_S=(I-(1/(double)dim)*I1*inv_C);

	tensor_S*= std::pow(J,(-2/(double)dim));

	double factor=(2*A*std::exp(-B*(I1_-3))+2*C*(I_m-3)/(I_m-I1_));

	tensor_S*= factor;
}


template <int dim>
void
IncompressibleYeohFleming<dim>::stress_S_vol(SymmetricTensor<2,dim> &tensor_S,
		  	  	  	  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	 double J=1;

	 J=sqrt(determinant(tensor_C));

	Assert(J>0,ExcMessage("Determinant J is smaller or equal to 0"));

	tensor_S=J*(B_vol/alpha)*std::sinh(alpha*(J-1))*invert(tensor_C);
}


template <int dim>
void
IncompressibleYeohFleming<dim>::material_tangent_iso(SymmetricTensor<4,dim> &tangent,
							  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{

	double J;

	J=sqrt(determinant(tensor_C));

	Assert(J>0,ExcMessage("Determinant J is smaller or equal to 0"));

	SymmetricTensor<2,dim> inv_C= invert(tensor_C);

	int I1=first_invariant(tensor_C);

	int I1_=std::pow(J,(-2/(double)dim))*I1;

	SymmetricTensor<2,dim> I;

	I=unit_symmetric_tensor<dim>();

	SymmetricTensor<4,dim> C_dev;

	double factor_dev=4*std::pow(J,(-2/(double)dim))*(-A*B*std::exp(-B*(I1_-3))+C*(I_m-3)/std::pow(I_m-I1_,2));

	double factor=(2/(double)dim)*(2*A*std::exp(-B*(I1_-3))+2*C*(I_m-3)/(I_m-I1_));

	C_dev=(product(I,I)-(1/(double)dim)*I1*product(I,inv_C)-(1/(double)dim)*I1*product(inv_C,I)+std::pow((1/(double)dim)*I1,2)*product(inv_C,inv_C));

	C_dev*=factor_dev;

	tangent=(tensorproduct(inv_C,inv_C)-(1/(double)dim)*product(inv_C,inv_C));

	tangent-=(product(inv_C,I)+product(I,inv_C)-(2/(double)dim)*I1_*product(inv_C,inv_C));

	tangent*=factor;

	tangent+=C_dev;

	tangent*=std::pow(J,(-2/(double)dim));
}


template <int dim>
void
IncompressibleYeohFleming<dim>::material_tangent_vol(SymmetricTensor<4,dim> &tangent,
							  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	 double J=1;

	 J=sqrt(determinant(tensor_C));

	Assert(J>0,ExcMessage("Determinant J is smaller or equal to 0"));

	SymmetricTensor<4,dim> C_vol1;
	SymmetricTensor<4,dim> C_vol2;

	SymmetricTensor<2,dim> C_inv= invert(tensor_C);

	C_vol1=product(C_inv,C_inv);

	C_vol1*=((B_vol/alpha)*J*std::sinh(alpha*(J-1))+std::pow(J,2)*B_vol*std::cosh(alpha*(J-1)));

	C_vol2=tensorproduct(C_inv,C_inv);

	C_vol2*=2*J*(B_vol/alpha)*std::sinh(alpha*(J-1));

	tangent=C_vol1-C_vol2;
}


}//namespace ConstitutiveLaws


namespace boost
{
	namespace serialization
	{
		template<class Archive, int dim>
		inline void save_construct_data(Archive & ar,
										const ConstitutiveLaws::IncompressibleYeohFleming<dim> *object_YF,
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
										ConstitutiveLaws::IncompressibleYeohFleming<dim> *object_YF,
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
			::new(object_YF)ConstitutiveLaws::IncompressibleYeohFleming<dim>(A,B,C,I_m,alpha,B_vol);
		}
	}//namespace serialization
}//namespace boost


#include "IncompressibleYeohFleming.inst"
