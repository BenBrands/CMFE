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

#pragma once


#include <deal.II/base/parameter_handler.h>


using namespace dealii;


namespace FE
{

class FESystemParameters
{
public:

	unsigned int poly_degree;

	static void declare_parameters(ParameterHandler &prm);

	void parse_parameters(ParameterHandler &prm);
};


class GeometryParameters
{
public:

	unsigned int global_refinement;
	double scale;
	double p_0;
	double parabola_root;

	static void declare_parameters(ParameterHandler &prm);

	void parse_parameters(ParameterHandler &prm);
};


class LinearSolverParameters
{
public:

	std::string type_lin;
	double tol_lin;
	double max_iterations_lin;
	std::string preconditioner_type;
	double preconditioner_relaxation;

	static void declare_parameters(ParameterHandler &prm);

	void parse_parameters(ParameterHandler &prm);
};


class NonlinearSolverParameters
{
public:

	unsigned int max_iterations_NR;
	double tol_f;
	double tol_u;

	static void declare_parameters(ParameterHandler &prm);

	void parse_parameters(ParameterHandler &prm);
};


class LoadSteppingParameters
{
public:

	unsigned int max_n_load_steps;
	double initial_load_factor;

	static void declare_parameters(ParameterHandler &prm);

	void parse_parameters(ParameterHandler &prm);
};


class FEParameters : public FESystemParameters,
                     public GeometryParameters,
                     public LinearSolverParameters,
                     public NonlinearSolverParameters,
					 public LoadSteppingParameters

{
public:

	FEParameters(const std::string &input_file);

	static void declare_parameters(ParameterHandler &prm);

	void parse_parameters(ParameterHandler &prm);
};

}
