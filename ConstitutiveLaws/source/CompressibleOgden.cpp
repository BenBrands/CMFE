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
#include <array>



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



template<int dim>
void tensor_to_symmetrictensor(SymmetricTensor<4,dim> &Sym, Tensor <4,dim> &T)
{
	for (unsigned int i=0; i<dim; i++)
			{
		for (unsigned int j=0; j<dim; j++)
					{
			for (unsigned int k=0; k<dim; k++)
						{
				for (unsigned int l=0; l<dim; l++)
							{
								Sym[i][j][k][l]=T[i][j][k][l];
							}
						}
					}
				}



}



template <int dim>
CompressibleOgden<dim>::CompressibleOgden(const Vector<double> alpha_,
													const Vector<double> mu_,
													const double kappa_vol_,
													const double betha_vol_,
													const double alpha_vol_,
													const double B_)
:
alpha(alpha_),
mu(mu_),
kappa_vol(kappa_vol_),
betha_vol(betha_vol_),
alpha_vol(alpha_vol_),
B_vol(B_)

{}



template <int dim>
double
CompressibleOgden<dim>::stress_S_iso_i(const std::array<std::pair<double,Tensor<1,dim>>,dim> &eigenpairs_C,
		  	  	  	  	  	  	  	   const double det_F,
									   const int i) const
{
	double result=0;

	const double lambda_i_bar = std::sqrt(eigenpairs_C[i].first) * std::pow(det_F,-1.0/(double)dim);

	for (unsigned int k=0; k<alpha.size(); ++k)
		{
			result += mu[k] * std::pow(lambda_i_bar,(alpha[k]))/ (eigenpairs_C[i].first );
		}

	for (unsigned int j=0; j<dim; ++j)
		{
			for (unsigned int k=0; k<alpha.size(); ++k)
				{
				const double lambda_j_bar = std::sqrt(eigenpairs_C[j].first) * std::pow(det_F,-1.0/(double)dim);

				result -= mu[k] * std::pow(lambda_j_bar,alpha[k]) / (eigenpairs_C[i].first* dim);
				}
		}
	return result;
}



template <int dim>
double
CompressibleOgden<dim>::dS_iso_i_dlambda_j(const std::array<std::pair<double,Tensor<1,dim>>,dim> &eigenpairs_C,
		  	  	  	  	  	  	  	   	   const double det_F,
										   const int i,const int j) const
{
	const double lambda_i_bar = std::sqrt(eigenpairs_C[i].first) * std::pow(det_F,-1.0/(double)dim);
	const double lambda_j_bar = std::sqrt(eigenpairs_C[j].first) * std::pow(det_F,-1.0/(double)dim);

	double result=0;

	if (i==j)
	{
		double term1=0;
		double term2=0;

		for (unsigned int p=0; p<alpha.size(); ++p)
			{
				double sum_k=0;
				for (unsigned int k=0; k<dim; ++k)
					{
						const double lambda_k_bar = std::sqrt(eigenpairs_C[k].first) * std::pow(det_F,-1.0/(double)dim);

						sum_k+= std::pow(lambda_k_bar,alpha[p])*(2.0+(alpha[p]/dim));
					}
				term1+= mu[p]*std::pow(lambda_i_bar,alpha[p])*((1.0-(double)2.0/dim)*alpha[p]-2.0);
				term2+=mu[p]*sum_k;
			}

		result=(term1/(eigenpairs_C[i].first*eigenpairs_C[i].first))+(term2/((double)dim*eigenpairs_C[i].first*eigenpairs_C[i].first));
	}
	else
	{
		double term1=0;

		for (unsigned int p=0; p<alpha.size(); ++p)
			{
				double sum_k=0;
				for (unsigned int k=0; k<dim; ++k)
					{
						const double lambda_k_bar = std::sqrt(eigenpairs_C[k].first) * std::pow(det_F,-1.0/(double)dim);

						sum_k+=std::pow(lambda_k_bar,alpha[p]);
					}
				term1+=(mu[p]*alpha[p]*((sum_k/(dim*dim))-(std::pow(lambda_i_bar,alpha[p])/dim)-(std::pow(lambda_j_bar,alpha[p])/dim)));
			}
		result=(term1/(eigenpairs_C[i].first*eigenpairs_C[j].first));
	}
	return result;
}



template <int dim>
void
CompressibleOgden<dim>::stress_S(SymmetricTensor<2,dim> &tensor_S,
									  const SymmetricTensor<2,dim> &tensor_C) const
{
	tensor_S=0;

	SymmetricTensor<2,dim> S_iso;

	SymmetricTensor<2,dim> S_vol;

	stress_S_iso(S_iso,tensor_C);

	stress_S_vol(S_vol,tensor_C);

	tensor_S+=S_iso;

	tensor_S+=S_vol;

}



template <int dim>
void
CompressibleOgden<dim>::material_tangent(SymmetricTensor<4,dim> &tangent,
									  	  	  const SymmetricTensor<2,dim> &tensor_C) const
{
	tangent=0;

	SymmetricTensor<4,dim> tangent_iso;

	SymmetricTensor<4,dim> tangent_vol;

	material_tangent_iso(tangent_iso,tensor_C);

	material_tangent_vol(tangent_vol,tensor_C);

	tangent+=tangent_iso;

	tangent+=tangent_vol;

}



template <int dim>
void
CompressibleOgden<dim>::stress_S_iso(SymmetricTensor<2,dim> &tensor_S,
		  	  	  	  	  	  		 const SymmetricTensor<2,dim> &tensor_C) const
{
	tensor_S = 0;

	std::array<std::pair<double,Tensor<1,dim>>,dim> eigenpairs = eigenvectors(tensor_C,SymmetricTensorEigenvectorMethod::ql_implicit_shifts);

	double det_F = 1;
	for (unsigned int i=0; i<dim; ++i)
		det_F *= eigenpairs[i].first;

	det_F = sqrt(det_F);

	for(unsigned int i=0; i<dim; i++)
		{
			double S_iso_i = stress_S_iso_i(eigenpairs,det_F,i);

			SymmetricTensor<2,dim> temp_S_iso=symmetrize(outer_product(eigenpairs[i].second,eigenpairs[i].second));

			temp_S_iso*=S_iso_i;

			tensor_S+=temp_S_iso;
		}
}



template <int dim>
void
CompressibleOgden<dim>::stress_S_vol(SymmetricTensor<2,dim> &tensor_S,
		  	  	  	  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	tensor_S=0;

	double J= sqrt(determinant(tensor_C));

	tensor_S = invert(tensor_C);

	tensor_S *= J*(B_vol/alpha_vol)*std::sinh(alpha_vol*(J-1));
}


template <int dim>
void
CompressibleOgden<dim>::material_tangent_iso(SymmetricTensor<4,dim> &tangent,
							  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	tangent=0;

	std::array<std::pair< double, Tensor< 1, dim>>,dim> eigenpairs= eigenvectors(tensor_C,SymmetricTensorEigenvectorMethod::ql_implicit_shifts);

	double det_F = 1;
	for (unsigned int i=0; i<dim; ++i)
		det_F *= eigenpairs[i].first;

	det_F = sqrt(det_F);

	Tensor<4,dim> temp_S_dev;


	for(unsigned int i=0; i<dim; ++i)
		{
		for(unsigned int m=0; m<dim; ++m)
			{
			Tensor<2,dim> product_ii=outer_product(eigenpairs[i].second,eigenpairs[i].second);

			Tensor<2,dim> product_mm=outer_product(eigenpairs[m].second,eigenpairs[m].second);

			double dS_iso_i_dlambda_m= dS_iso_i_dlambda_j(eigenpairs, det_F, i, m);

			temp_S_dev+=(outer_product(product_ii,product_mm)*dS_iso_i_dlambda_m);
			}
		}

	Tensor<4,dim> S;
	Tensor<4,dim> S_i_j_diff;
	Tensor<4,dim> S_i_j_same;



	for(unsigned int l=0;l<dim;l++)
		{
			for(unsigned int k=0;k<dim;k++)
				{
					Tensor<2,dim> product_lk=outer_product(eigenpairs[l].second,eigenpairs[k].second);

					Tensor<2,dim> product_kl=outer_product(eigenpairs[k].second,eigenpairs[l].second);

					double temp_S_i_j_diff=0;

					double temp_S_i_j_same=0;

					if (l!=k)
					{
						if(std::fabs(std::sqrt(eigenpairs[k].first)-std::sqrt(eigenpairs[l].first))>1e-8)
						{
						double S_k=stress_S_iso_i(eigenpairs,det_F,k);

						double S_l=stress_S_iso_i(eigenpairs,det_F,l);

						temp_S_i_j_diff=(S_k-S_l)/(eigenpairs[k].first-eigenpairs[l].first);
						}
						else
						{
							double dSk_dlambdak=dS_iso_i_dlambda_j(eigenpairs,det_F,k,k);

							double dSl_dlambdak=dS_iso_i_dlambda_j(eigenpairs,det_F,l,k);

							temp_S_i_j_same=(dSk_dlambdak-dSl_dlambdak)*0.5;
						}
					}
					S_i_j_diff+=outer_product(product_lk,product_lk)*temp_S_i_j_diff;

					S_i_j_diff+=outer_product(product_lk,product_kl)*temp_S_i_j_diff;

					S_i_j_same+=outer_product(product_lk,product_lk)*temp_S_i_j_same;

					S_i_j_same+=outer_product(product_lk,product_kl)*temp_S_i_j_same;
				}
			}

	S=S_i_j_diff+S_i_j_same;

	S+=temp_S_dev;

	SymmetricTensor<4,dim> S_4;

	tensor_to_symmetrictensor(S_4,S);

	tangent=S_4;

	S=0;

	S_4=0;
}



template <int dim>
void
CompressibleOgden<dim>::material_tangent_vol(SymmetricTensor<4,dim> &tangent,
							  	  	  			  const SymmetricTensor<2,dim> &tensor_C) const
{
	tangent=0;

	double J = sqrt(determinant(tensor_C));

	const double sinh_a_J = std::sinh(alpha_vol*(J-1));

	const double cosh_a_J = std::cosh(alpha_vol*(J-1));

	SymmetricTensor<4,dim> C_vol1;
	SymmetricTensor<4,dim> C_vol2;

	SymmetricTensor<2,dim> C_inv= invert(tensor_C);

	C_vol1=outer_product(C_inv,C_inv);

	double factor = (B_vol/alpha_vol)*J*sinh_a_J;

	factor += std::pow(J,2)*B_vol*cosh_a_J;

	C_vol1 *= factor;

	C_vol2=tensorproduct(C_inv,C_inv);

	C_vol2*=2*J*(B_vol/alpha_vol)*sinh_a_J;

	tangent=C_vol1-C_vol2;
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
			ar << object_O->kappa_vol;
			ar << object_O->betha_vol;
			ar << object_O->alpha_vol;
			ar << object_O->B_vol;
		}


		template<class Archive, int dim>
		inline void load_construct_data(Archive & ar,
										ConstitutiveLaws::CompressibleOgden<dim> *object_O,
										const unsigned int /*version*/)
		{
			// retrieve data from archive required to construct new instance
			Vector<double> alpha, mu;
			double alpha_vol=0, B_vol=0, kappa_vol=0, betha_vol=0;
			ar >> alpha;
			ar >> mu;
			ar >> kappa_vol;
			ar >> betha_vol;
			ar >> alpha_vol;
			ar >> B_vol;


			// invoke inplace constructor to initialize instance of class
			::new(object_O)ConstitutiveLaws::CompressibleOgden<dim>(alpha,mu,kappa_vol,betha_vol,alpha_vol,B_vol);
		}
	}//namespace serialization
}//namespace boost


#include "CompressibleOgden.inst"
