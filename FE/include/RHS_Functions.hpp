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


#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/lac/vector.h>


using namespace dealii;


namespace FE
{

template <int dim>
class Load : public Function<dim>
{
public:

	Load(const double lf=1);

	virtual ~Load() = default;

	virtual void set_load_factor(const double lf) {load_factor_=lf;}

	virtual double load_factor() const {return load_factor_;}

	virtual void tensor_value(const Point<dim> &p,
							  Tensor<1,dim> &t) const=0;

	virtual void tensor_value_list(const std::vector<Point<dim>> &point_list,
								   std::vector<Tensor<1,dim>> &value_list) const=0;

private:

	double load_factor_;
};


// formula: f(X;X_0) = - p_0 / (d/r)^2 + p_0
// d = || X - X_0 || is distance of X from origin X_0; r is the distance where the load changes its sign
// force is acting in negative direction of last coordinate axis: y-direction in 2D and z-direction in 3D
template <int dim>
class ParabolicLoad : public Load<dim>
{
public:

	ParabolicLoad(const Point<dim> &origin,
				  const double p_0,
				  const double radius,
				  const double lf=1);

	virtual ~ParabolicLoad() = default;

	virtual double value(const Point<dim> &p,
						 const unsigned int component=0) const;

	virtual void vector_value(const Point<dim> &p,
							  Vector<double> &value) const;

	virtual void value_list(const std::vector<Point<dim>> &point_list,
							std::vector<double> &value_list,
							const unsigned int component=0) const;

	virtual void vector_value_list(const std::vector<Point<dim>> &point_list,
								   std::vector<Vector<double>> &value_list) const;

	virtual void tensor_value(const Point<dim> &p,
							  Tensor<1,dim> &t) const;

	virtual void tensor_value_list(const std::vector<Point<dim>> &point_list,
								   std::vector<Tensor<1,dim>> &value_list) const;

private:

	const Point<dim> origin;

	const double p_0;

	const double radius;
};


}//namespace FE
