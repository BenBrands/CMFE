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


#include "CompressibleOgden.hpp"


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

template<typename T>
int kroneckerdelta(T &x,T &y)
{
	if (x==y)
		return 1;
	else
		return 0;

}

template <int dim>
CompressibleOgden<dim>::CompressibleOgden(const Vector<double> alpha_,
													const Vector<double> mu_)
:
alpha(alpha_),
mu(mu_)
{}

template <int dim>
void
CompressibleOgden<dim>::stress_S_i(double &S,
		  const SymmetricTensor<2,dim> &tensor_C,const int i) const
{
double J=1;

J=sqrt(determinant(tensor_C));

Assert(J>0,ExcMessage("Determinant J is smaller or equal to 0"));

double S_1;

double S_2;

Vector<double> lambda_(dim);

for (unsigned int c=0; c<dim; c++)
	 {
		 lambda_[c]=sqrt(tensor_C[c][c]);
	 }


lambda_*=std::pow(J,-1/(double)dim);

S_1=0;

S_2=0;

	for (unsigned int k=0; k<alpha.size(); k++)
	{
		S_1+=std::pow(lambda_[i],alpha[k]-2);

		S_1*=mu[k];



		for (unsigned int j=0;j<dim; j++)
		{
			S_2+=(mu[k]/(double)dim);

			S_2*=std::pow(lambda_[j],alpha[k]);

			S_2*=std::pow(lambda_[i],-2);

		}
	}

S=S_1;

S-=S_2;

S*=std::pow(J,-2/(double)dim);
}

template <int dim>
void
CompressibleOgden<dim>::stress_S(SymmetricTensor<2,dim> &tensor_S,
									  const SymmetricTensor<2,dim> &tensor_C) const
{

	for(unsigned int i=0; i<dim; i++)
	{
		double S=0;

		stress_S_i(S,tensor_C,i);

		tensor_S[i][i]=S;

	}

}


template <int dim>
void
CompressibleOgden<dim>::material_tangent(SymmetricTensor<4,dim> &tangent,
									  	  	  const SymmetricTensor<2,dim> &tensor_C) const
{

	 double J=1;

	 J=sqrt(determinant(tensor_C));

	Assert(J>0,ExcMessage("Determinant J is smaller or equal to 0"));

	Vector<double> lambda_(dim);

	 for (unsigned int c=0; c<dim; c++)
		 {
			 lambda_[c]=sqrt(tensor_C[c][c]);
		 }

	lambda_*= std::pow(J,-1/(double)dim);

	SymmetricTensor<4,dim> dev_S;

	SymmetricTensor<4,dim> S;


	for(unsigned int i=0; i<dim; i++)
	{
		for(unsigned int m=0; m<dim; m++)
			{

			dev_S[i][i][m][m]=0;

			for(unsigned int k=0; k<alpha.size(); k++)
				{
					int kroneckerdeltaIM= kroneckerdelta(i,m);

					dev_S[i][i][m][m]+=mu[k]*((-2/(double)dim)*std::pow(lambda_[m],-2)*std::pow(lambda_[i],alpha[k]-2)+(alpha[k]-2)*(kroneckerdeltaIM*std::pow(lambda_[m],-1)*std::pow(lambda_[i],alpha[k]-3)-(1/(double)dim)*std::pow(lambda_[i],alpha[k]-2)*std::pow(lambda_[m],-2))-(alpha[k]/(double)dim)*std::pow(lambda_[i],-2)*std::pow(lambda_[m],alpha[k]-2));

					for (unsigned int j=0; j<dim; j++)
					{
						dev_S[i][i][m][m]+=(mu[k]/(double)dim)*((2/(double)dim)*std::pow(lambda_[j],alpha[k])*std::pow(lambda_[i],-2)*std::pow(lambda_[m],-2)-alpha[k]*(-(1/(double)dim)*std::pow(lambda_[j],alpha[k])*std::pow(lambda_[i],-2)*std::pow(lambda_[m],-2))+2*(kroneckerdeltaIM*std::pow(lambda_[j],alpha[k])*std::pow(lambda_[i],-3)*std::pow(lambda_[m],-1)-(1/(double)dim)*std::pow(lambda_[i],-2)*std::pow(lambda_[j],alpha[k])*std::pow(lambda_[m],-2)));
					}
				}

			}
	}

			for(unsigned int l=0;l<dim;l++)
			{
				for(unsigned int k=0;k<dim;k++)
				{

					if ((l!=k)&& (lambda_[l]!=lambda_[k]))
					{
					double factor;

					factor=std::pow(lambda_[k],2);

					factor-=std::pow(lambda_[l],2);

					double S_k=0;

					double S_l=0;

					stress_S_i(S_k,tensor_C,k);

					S[l][k][l][k]=S_k;

					S[l][k][k][l]=S_k;

					stress_S_i(S_l,tensor_C,l);

					S[l][k][l][k]-=S_l;

					S[l][k][k][l]-=S_l;

					tangent[l][k][l][k]+=S[l][k][l][k];

					tangent[l][k][k][l]+=S[l][k][k][l];

					tangent[l][k][l][k]/=factor;

					tangent[l][k][k][l]/=factor;

}
					else
					{
						tangent[l][k][l][k]=0;

						tangent[l][k][k][l]=0;
					}

				}
			}


		dev_S*=std::pow(J,-(4/(double)dim));

		tangent+=dev_S;

}


template <int dim>
void
CompressibleOgden<dim>::stress_S_iso(SymmetricTensor<2,dim> &tensor_S,
		  	  	  	  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	(void)tensor_S;
	(void)tensor_C;

	Assert(false,ExcMessage("The function stress_S_iso is not callable for CompressibleOgden"));
}


template <int dim>
void
CompressibleOgden<dim>::stress_S_vol(SymmetricTensor<2,dim> &tensor_S,
		  	  	  	  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	(void)tensor_S;
	(void)tensor_C;

	Assert(false,ExcMessage("The function stress_S_vol is not callable for CompressibleOgden"));
}


template <int dim>
void
CompressibleOgden<dim>::material_tangent_iso(SymmetricTensor<4,dim> &tangent,
							  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	(void)tangent;
	(void)tensor_C;

	Assert(false,ExcMessage("The function material_tangent_iso is not callable for CompressibleOgden"));
}


template <int dim>
void
CompressibleOgden<dim>::material_tangent_vol(SymmetricTensor<4,dim> &tangent,
							  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	(void)tangent;
	(void)tensor_C;

	Assert(false,ExcMessage("The function material_tangent_vol is not callable for CompressibleOgden"));
}


}//namespace ConstitutiveLaws


namespace boost
{
	namespace serialization
	{
		template<class Archive, int dim>
		inline void save_construct_data(Archive & ar,
										const ConstitutiveLaws::CompressibleOgden<dim> *object_O,
										const unsigned int /*version*/)
		{
			// save data required to construct instance
			ar << object_O->alpha;
			ar << object_O->mu;
		}


		template<class Archive, int dim>
		inline void load_construct_data(Archive & ar,
										ConstitutiveLaws::CompressibleOgden<dim> *object_O,
										const unsigned int /*version*/)
		{
			// retrieve data from archive required to construct new instance
			Vector<double> alpha, mu;
			ar >> alpha;
			ar >> mu;

			// invoke inplace constructor to initialize instance of class
			::new(object_O)ConstitutiveLaws::CompressibleOgden<dim>(alpha,mu);
		}
	}//namespace serialization
}//namespace boost


#include "CompressibleOgden.inst"
