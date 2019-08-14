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


#include <iostream>
#include <fstream>


#include "SolidSerial.hpp"
#include "RHS_Functions.hpp"


int main(int argc, char **argv)
{
	using namespace dealii;

	try
	{
		AssertThrow(argc>2,ExcMessage("Two parameter file names have to be passed!"));

		std::vector<std::string> args;

		for (int i=1; i<argc; ++i)
			args.push_back (argv[i]);

		const std::string fe_parameter_file = args[0];

		const std::string material_parameter_file = args[1];

		const unsigned int dim = 3;

		Triangulation<dim> triangulation;

		FE::SolidSerial<dim> block(triangulation,
								   fe_parameter_file,
								   material_parameter_file);

		block.run();
	}
	catch (std::exception &exc)
    {
		std::cerr << std::endl
				  << std::endl
				  << "----------------------------------------------------"
				  << std::endl;
		std::cerr << "Exception on processing: " << std::endl
                  << exc.what() << std::endl
				  << "Aborting!" << std::endl
				  << "----------------------------------------------------"
				  << std::endl;

		return 1;
    }
	catch (...)
	{
		std::cerr << std::endl
				  << std::endl
				  << "----------------------------------------------------"
				  << std::endl;
		std::cerr << "Unknown exception!" << std::endl
				  << "Aborting!" << std::endl
				  << "----------------------------------------------------"
				  << std::endl;
		return 1;
	}
	return 0;
}
