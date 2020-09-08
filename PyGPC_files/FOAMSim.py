#!/usr/bin/python3
# -*- coding: utf-8 -*-
import re
import os
import numpy as np
import csv

from subprocess import check_output
from subprocess import run
from shutil import copytree
from shutil import rmtree
from shutil import move

from pygpc.AbstractModel import AbstractModel

'''
 This class extents the "Model" abstract class of the PyGPC framework.
 It encapsulates the handling of the OpenFOAM test case:
    - setting the conductivity values of this particular simualation pass
    - initiating the simulation 
    - retrieving the simulation results
 
 Required parameters = OF_CASE_DIR_BASE :
    - The base directory of a readily setup OpenFOAM case.
    - The case directory must already contain a valid mesh and
      boundary conditions.

 @author: Benjamin Kalloch
'''
class FOAMSim(AbstractModel):

    def simulate(self, process_id, matlab_engine):
        of_case_dir = self.of_cas_dir_base + "_" + str(process_id)

        # Step 1: setup the OpenFOAM simulation case by replacing the conductivities with
        #         the values provided ny the PyGPC framework
        rmtree(of_case_dir + "/0", ignore_errors=True)
        copytree(of_case_dir + "/0_clean", of_case_dir + "/0")

        with open(of_case_dir + "/0/sigma") as f:
            sigma = f.read()

        sigma = sigma.replace("%SKIN_VAL%", str(self.p["scalp_cond"][0]))
        sigma = sigma.replace("%SKULL_VAL%", str(self.p["skull_cond"][0]))
        sigma = sigma.replace("%CSF_VAL%", str(self.p["csf_cond"][0]))
        sigma = sigma.replace("%GM_VAL%", str(self.p["gm_cond"][0]))
        sigma = sigma.replace("%WM_VAL%", str(self.p["wm_cond"][0]))
        sigma = sigma.replace("%AIR_VAL%", str(self.p["air_cond"][0]))
        sigma = sigma.replace("%ELECTRODE_VAL%", str(self.p["electrode_cond"][0]))
        
        if "lesion_cond" in self.p:
            sigma = sigma.replace("%LESION_VAL%", str(self.p["lesion_cond"][0]))

        with open(of_case_dir + "/0/sigma", "w") as f:
            f.write(sigma)

        # Step 2: Run the simulation
        #         The simulation is execute using a bash script, which is located in the base dir.
        #         We provide its location as an argument to pass on to the solver application.
        stdout_string = check_output(["bash " + of_case_dir + "/runSim.sh " + of_case_dir], shell=True)
        with open(of_case_dir + "/solver.log", "a") as f:
            f.write("\nSolver output of iteration #"+str(self.i_grid)+":\n")
            f.write(str(stdout_string.decode("utf-8")))

            # we must check the number of iterations, because the simulation results will be stored in a
            # directory named numerically "#iterations + 1"
        regex_result   = re.search("SIMPLE solution converged in ([0-9]+) iterations", stdout_string.decode("utf-8"))
        num_iterations = int(regex_result.group(1))

        # Step 3: Query the results
        internal_field_flag = False
        num_lines_to_read   = -1
        ElPot = np.array([])

        result_dir = of_case_dir + "/" + str( num_iterations + 1);

        with open(result_dir + "/ElPot") as f:
            for line in f:
                if "internalField" in line:
                    num_lines_to_read = int(next(f))
                    next(f)
                    break
            ElPot = np.zeros(num_lines_to_read, dtype="float64")

            for i in range(0, num_lines_to_read):
                ElPot[i] = float(next(f))

        # sample the electrical field onto the midlayer_surface
        script_location=of_case_dir + "/interpolate_to_midlayer_and_retrieve_values.sh"
        param_surface=of_case_dir+"/midlayer.stl"
        param_casefile=of_case_dir + "/dummy.OpenFOAM"
        param_outfile=result_dir+"/E_midlayer.csv"
        # create dummy OpenFOAM case-file for paraview
        open(param_casefile, 'a').close()

        check_output(["bash " + script_location + " " + param_surface + " " + param_casefile + " " + param_outfile], shell=True)
        '''
        ############################################################
        # validation occurs for E-midlayer only
        # REMOVE THIS FOR NORMAL SENSITIVITY ANALYSIS!!!
        # This is only needed for MC validation for the GPC
        midlayer_data = []
        with open( of_case_dir + "/result_" + str(process_id) + "_iter_6/E_midlayer.csv" ) as midlayer_data_csv:
            reader = csv.reader( midlayer_data_csv, delimiter=',')
            next(reader) # skip first line
            for row in reader:
                midlayer_data.append( float(row[0]) )
        ###########################################################
        '''
        # Optional: rename the result directory, appending the iteration number to save the content for later analysis
        move(result_dir, of_case_dir + "/result_" + str(self.i_grid) + "_iter_" + str(num_iterations + 1))

        # OpenFOAM also creates a directory "#iterations". We do not need this directory
        # but it may confuse ParaView concerning the timesteps if different runs need different
        # amounts of iterations. Thus, we remove it.
        rmtree(of_case_dir + "/" + str(num_iterations), ignore_errors=True)
        
        return ElPot

    # Do nothing since there is no combination of input paramters the simulation cannot perform.
    def validate(self):
        pass

    def set_parameters( self, p, context ):
        self.of_cas_dir_base = p["OF_CASE_DIR_BASE"]
        return super().set_parameters(p, context)

    '''
     This static method is not part of the abstract model of the PyGPD framework.
     It is used to write the PyGPC results (variance, Sobol indices) in the OpenFOAM
     file format. 
    '''
    @staticmethod
    def write_result_field(case_dir, fieldName, data):
        out_dir = case_dir + "/999/"

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        with open(case_dir + "/Field_template") as f:
            template = f.read()

        data_str = np.char.mod('%f', data)
        data_str = "\n".join(data_str)

        template = template.replace("%FIELDNAME%", fieldName)
        template = template.replace("%NUM_VALUES%", repr(len(data)))
        template = template.replace("%DATA%", data_str)

        with open(out_dir + "/" + fieldName, "w") as f:
            f.write(template)
