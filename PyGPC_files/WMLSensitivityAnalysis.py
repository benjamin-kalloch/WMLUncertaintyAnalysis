#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
This script implements the PyGPC framework to perform a senstivity analysis
on tDCS simulations including information of white matter lesion tissue.

@author: Benjamin Kalloch
"""
import pygpc
import h5py
import numpy as np
import csv
import sys, getopt

import pickle

from benjamin.FOAMSim import FOAMSim

from re import search
from itertools import chain
from collections import OrderedDict
from pygpc.Problem import Problem
from pygpc.RandomParameter import Beta
from os import walk

class WMLSensitivityAnalysis:
    def __init__(self, case_dir_base, results_filename):
        self.results_filename = results_filename
        self.of_case_dir_base = case_dir_base

    def run(self, withLesion=False):
        # define the Problem
        parameters = OrderedDict()  # conductivity values adapted according to (Saturnino et al. 2019)
        parameters["scalp_cond"]     = Beta(pdf_shape=[3, 3], pdf_limits=[.2, .5])     # before: [.28,.87]
        parameters["skull_cond"]     = Beta(pdf_shape=[3, 3], pdf_limits=[.003, .012]) # compact bone, before: [.0016,.033] 
        parameters["csf_cond"]       = Beta(pdf_shape=[3, 3], pdf_limits=[1.2, 1.8])   # before: 1.65 (fixed not random variable)
        parameters["gm_cond"]        = Beta(pdf_shape=[3, 3], pdf_limits=[.1, .6])     # before: [.22, .67]
        parameters["wm_cond"]        = Beta(pdf_shape=[3, 3], pdf_limits=[.1, .4])     # before: [.09, .29]
        parameters["air_cond"]       = [1e-15]
        parameters["electrode_cond"] = [1.4]
        if withLesion:
            print("Performing sensitivity analysis with lesioned white matter tissue as input variable.")
            parameters["lesion_cond"] = Beta(pdf_shape=[3, 3], pdf_limits=[.1, 1.8])    # range: [lowest-WM,highest-CSF], before: [.04, 1.5]
        else:
            print("Performing sensitivity analysis without lesioned white matter tissue.")

        parameters["OF_CASE_DIR_BASE"] = self.of_case_dir_base

        wmlSensitivity = Problem(FOAMSim(), parameters)

        # gPC options
        options = dict()
        options["order_start"]       = 1     # default
        options["matrix_ratio"]      = 1.5   # default
        options["eps"]               = 1e-3  # default 
        options["order_max_norm"]    = 1     # default
        options["order_end"]         = 5
        options["interaction_order"] = 3     # Auf 2 oder 3 setzen
        options["seed"]              = 42
        options["n_cpu"]             = 6
        options["settings"] = None  # ?
        options["fn_results"] = self.results_filename
        options["solver"]     = "Moore-Penrose" # default
        options["error_norm"] = "relative"      # default
        options["method"]     = "Reg"           # not necessary when using "RegAdaptive" ?

        options["print_func_time"]   = True     # default
        options["gradient_enhanced"] = False    # default
        options["projection"]        = False    # default
        options["adaptive_sampling"] = False    # default
        options["GPU"]               = False    # default
        
        # define the kind of gPC algorithm
        algorithm = pygpc.RegAdaptive(problem=wmlSensitivity, options=options)

        # run gPC algorithm
        print( "(1) RUN GPC" )
        gpc_phi, coeffs_phi, results_phi = algorithm.run()
        ######### DEBUG START
        '''
        print("Loading:", self.results_filename+".hdf5")
        print("Loading:", self.results_filename+".pkl")
        results_phi = []
        try:
            with h5py.File( self.results_filename+".hdf5","r") as f:
               results_phi = f['model_evaluations/results'][:]
        except (KeyError, ValueError):
            return None

        gpc_phi=None
        with open( self.results_filename+".pkl", "rb") as p:
            gpc_phi=pickle.load(p)
        '''
        ######### DEBUG END    
        # post process
        print( "(2) COMPUTE SOBOL INDICES" )
        pygpc.get_sensitivities_hdf5(fn_gpc=options["fn_results"],
                                     output_idx=None,       # Indices of output quantities (QOIs) to consider in postprocessind
                                     calc_sobol=True,       # Calculate Sobol indices
                                     calc_global_sens=True, # Calculate global derivative based sensitivities
                                     calc_pdf=False)         # Calculate probability density functions of output quantities
        print( "(3) WRITE POSTPORCESSED DATA AS OF FIELD" )
        structures = np.array(["scalp", "skull", "csf", "gm", "wm"])
        if withLesion:
            structures.append("lesion")

        print("Output will contain Sobol Indices of the following structures:", structures)

        # read hdf5 files with post-processed results and write as OpenFoam field
        with h5py.File(self.results_filename + ".hdf5", "r") as f:
            try:
                FOAMSim.write_result_field(parameters["OF_CASE_DIR_BASE"], "ElPot_mean", f["sens/mean"][0])
                FOAMSim.write_result_field(parameters["OF_CASE_DIR_BASE"], "ElPot_var", f["sens/var"][0])
                FOAMSim.write_result_field(parameters["OF_CASE_DIR_BASE"], "ElPot_rstd", f["sens/rstd"][0])
                # "sobol_idx_order" is a boolean array representing which of the random input variables
                # is represented by the data on "sens/sobol" of the same index.
                # e.g. given 3 input random variables A,B,C, sobol_idx_bool = [True, False, False]
                #      means that sobol[0] represents the contribution of A to the variance
                for sobol_order in range(0, len( f["sens/sobol_idx_bool"] ) ):
                    sobol_structures = structures[ np.array( f["sens/sobol_idx_bool"][sobol_order] ) ]
                    structures_string = ""
                    for structure in sobol_structures:
                        structures_string += structure + "_"
                    structures_string = structures_string[:-1]
                    FOAMSim.write_result_field(parameters["OF_CASE_DIR_BASE"], "ElPot_sobol_"+str(sobol_order)+"_"+structures_string, f["sens/sobol"][sobol_order])
            except KeyError:
                print("Cannot read result data.")

        # Determine gPC coefficients again for the E field on the midlayer.
        #   -> We use the before determined expansion (from ElPot).
        #   -> We determine the coefficients based on the E field sampled on the midlayer.
        total_num_iterations = len( results_phi )
        #total_num_iterations = 51

        E_midlayer = [0]*total_num_iterations
        
        # The results are each written into a separate directory.
        # Each process creates one sub-case-folder and stores its results in that folder.
        # The results of each process are stored in a sub-directory of the corresponding case folder.
        # This result directory is names results_'Iterationnumber'_iter_'NumderOfOpenFOAMIterationsToConverge'
        #   e.g. result_0_iter8 = is the result of the very first pygpc-iteration and OpenFOAM took 8
        #        iterations to converge below the global residual.

        # 1st: collect the sub-process-case-directories
        process_directories=[ parameters["OF_CASE_DIR_BASE"]+"_0" ]
        for i in range(1,options["n_cpu"]):
            process_directories.append( parameters["OF_CASE_DIR_BASE"]+"_"+str(i))

        print ("(4) COLLECT RESULTS OF MIDLAYER" )
        # 2nd: collect the results from the respective result-folders of each sub-process-directory
        #      and add them to the result-array (E_midlayer) to the correct position (according to
        #      the iteration number when that result was calculated)
        for root, dirs, files in chain.from_iterable( walk(path) for path in process_directories ):
            iteration_number = search("result_([0-9]+)_iter.*", root)
            if iteration_number:
                with open( root + "/E_midlayer.csv") as midlayer_data_csv:
                    midayer_data = []
                    reader = csv.reader( midlayer_data_csv, delimiter=',')
                    next(reader) # skip first line
                    for row in reader:
                        midayer_data.append( float(row[0]) )
                    E_midlayer[ int(iteration_number.group(1)) ] = midayer_data
        '''
        # DEBUG: for debugging the later steps without having to re-compute the gPC for Phi
        with open( self.results_filename+".pkl","rb") as pickle_file:
            gpc_phi = pickle.load( pickle_file ) 
        gpc_phi.n_cpu=6

        '''
        print( "(5) COMPUTE NEW GPC COEFFICIENTS FROM MIDLAYER DATA")
        # Compute new coefficients according to the E field at the midlayer
        # (+ store for later/repeated use)
        print( len( E_midlayer ) )
        print( len( E_midlayer[0] ) )
        coeffs_E = gpc_phi.solve(E_midlayer, solver=options["solver"])
     
     
        print( "(6) WRITE E MIDLAYER")
        print("Writing to '" + self.results_filename+"_coeffs_E.hdf'")   
        with h5py.File(self.results_filename+"_coeffs_E.hdf5", "a") as f:
            if "coeffs_E" in f.keys():
                del f['coeffs_E']
            f.create_dataset("coeffs_E", data=coeffs_E, maxshape=None, dtype="float64")
        '''

        # DEBUG: for debugging the later stages without having to re-calculate the above again...
        coeffs_E = []
        try:
            with h5py.File( self.results_filename+"_coeffs_E.hdf5","r") as f:
               coeffs_E = f['coeffs_E'][:]
        except (KeyError, ValueError):
            return None

        # determine error
        '''
        print( "(7) DETERMINE LOOCV ERROR")
        eps = gpc_phi.loocv(results=np.array(E_midlayer), coeffs=coeffs_E )   # I do not need a nonNanMask since my results won't become NaN
        '''
        pygpc.validate_gpc_mc( gpc_phi, coeffs_E, n_samples=int(6), n_cpu=6, fn_out=self.results_filename+"_MC_validation_result")
        exit()
        '''
        print("*** Error in E field estimation = " + str(eps) + " ***" );
        
        # postprocessing of E-field gPC
        mean_E  = gpc_phi.get_mean(coeffs_E[:, np.arange(coeffs_E.shape[1])])[0]    # in case of multi-part gpce this might be > 0
        
        std_E = gpc_phi.get_std(coeffs_E[:, np.arange(coeffs_E.shape[1])])[0]
        var_E = std_E ** 2
        
        sobol_E, sobol_idx_E, sobol_idx_bool_E = gpc_phi.get_sobol_indices( coeffs_E[:, np.arange(coeffs_E.shape[1]) ])
        structures_strings=[]
        for sobol_order in range(0, len( sobol_idx_bool_E ) ):
            sobol_structures = structures[ np.array( sobol_idx_bool_E[sobol_order] ) ]
            structures_string = ""
            for structure in sobol_structures:
                structures_string += structure + "_"
            structures_strings.append(structures_string[:-1])

        with open( self.results_filename+"_postprocessing_E.csv", "w", newline="" ) as midlayer_gpc_results_csv:
            writer = csv.writer( midlayer_gpc_results_csv, delimiter="," )
            sobol_indices_columnnames = []
            for sobol_order in range(0,len(sobol_E)):
                sobol_indices_columnnames.append( "E_sobol_"+str(sobol_order)+"_"+structures_strings[sobol_order] )

            writer.writerow(["E_mean", "E_var"] + sobol_indices_columnnames)
            for i in range(0, len(mean_E) ):
                sobol_indices = []
                for sobol_order in range(0, len(sobol_E)):
                    sobol_indices.append( sobol_E[sobol_order][i] )
                writer.writerow( [mean_E[i], var_E[i]] + sobol_indices )

def printHelpText():
    print( "WMLSensitivityAnalysis.py -b|--base_dir PATH_TO_SIM_BASE_FOLDER -o|--out_file PATH_TO_OUTPUT_DIRECTORY" )

def main(argv):
    try:
        opts, args = getopt. getopt( argv, "hb:o:",["base_dir=","out_file="])
    except getopt.GetoptError:
        printHelpText()
        sys.exit(2)

    simulation_base_folder = None
    output_file = None
    for opt, arg in opts:
        if opt == '-h':
            printHelpText()
            sys.exit()
        elif opt in ("-b", "--base_dir"):
            simulation_base_folder = arg
        elif opt in ("-o", "--out_file"):
            output_file = arg
            
    print("Simulation base directory:", simulation_base_folder)
    print("Output files:", output_file)
    if simulation_base_folder is not None and output_file is not None:
        sim = WMLSensitivityAnalysis(simulation_base_folder,output_file)
        sim.run(withLesion=False)
    else:
        print("Invalid input parameters!")
        printHelpText()

if __name__ == "__main__":
    main(sys.argv[1:])
