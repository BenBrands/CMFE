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


#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_point_data.h>
#include <deal.II/base/smartpointer.h>

#include <deal.II/grid/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/dofs/dof_handler.h>


#include <iostream>
#include <fstream>
#include <algorithm>


#include "SimulationParameters.hpp"
#include "QPData.hpp"
#include "RHS_Functions.hpp"


using namespace dealii;


namespace FE
{

template <int dim>
class SolidSerial
{
public:

	SolidSerial(Triangulation<dim> &tria,
				const std::string &fe_parameter_file,
				const std::string &material_parameter_file);

	void run();

private:

	void make_grid();

	void system_setup();

	void make_constraints(AffineConstraints<double> &constraints,
						  const unsigned int it_nr=0) const;

	void setup_qph(CellDataStorage<typename DoFHandler<dim>::cell_iterator,
				   QPData<dim>> &qp_data) const;

	bool update_qph(CellDataStorage<typename DoFHandler<dim>::cell_iterator,
					QPData<dim>> &qp_data,
					const Vector<double> &solution) const;

	void assemble_matrix(SparseMatrix<double> &matrix,
						 CellDataStorage<typename DoFHandler<dim>::cell_iterator,
						 	 	 	 	 QPData<dim>> &qp_data,
						 const AffineConstraints<double> &constraints) const;

	void assemble_residuum(Vector<double> &residuum,
			 	 	 	   CellDataStorage<typename DoFHandler<dim>::cell_iterator,
			 	 	 	 	 	   	   	   QPData<dim>> &qp_data,
						   std::shared_ptr<Load<dim>> load,
						   const AffineConstraints<double> &constraints) const;

	bool solve_NewtonRaphson(Vector<double> &solution,
							 CellDataStorage<typename DoFHandler<dim>::cell_iterator,
					 	 	 	 	 	   	 QPData<dim>> &qp_data,
							 std::shared_ptr<Load<dim>> load) const;

	std::pair<unsigned int,double>
	solve_linear_system(Vector<double> &solution,
						const SparseMatrix<double> &matrix,
						const Vector<double> &rhs) const;

	void output_results(const Vector<double> &vec,
						CellDataStorage<typename DoFHandler<dim>::cell_iterator, QPData<dim>> &qp_data,
						const std::string &file_name) const;

	void output_results(const Vector<double> &vec,
						const std::string &file_name) const;

	FEParameters parameters;

	SmartPointer<Triangulation<dim>> triangulation;

	const FESystem<dim> fe;

	DoFHandler<dim> dof_handler;

	const FEValuesExtractors::Vector u_fe;

	const FEValuesExtractors::Scalar u_x_fe, u_y_fe, u_z_fe;

	const QGauss<dim> quadrature;

	const QGauss<dim-1> face_quadrature;

	const MaterialData material;

	SparsityPattern sparsity_pattern;
};


}//namespace FE
