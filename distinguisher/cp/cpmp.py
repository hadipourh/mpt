#!/usr/bin/python3
#-*- coding: UTF-8 -*-
"""
Applying the monomial prediction technique to find integral
distinguishers of WARP in the single-key setting
Copyright (C) 2021  Hosein Hadipour

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os
import time
import minizinc
import datetime
from argparse import ArgumentParser, RawTextHelpFormatter

"""
x_roundNumber_nibbleNumber_bitNumber
x_roundNumber_nibbleNumber_0: msb
x_roundNumber_nibbleNumber_3: lsb
Variable mapping:

x_r_0                       k_r  x_r_1  ...
|                            |     |
|--- y_r---> | S | ---z_r----+---->+    ...
|                                  |    ...

"""

class Warp:
    count = 0
    def __init__(self, nrounds=10, cp_solver_name="or-tools", time_limit=-1):
        Warp.count += 1
        self.nrounds = nrounds
        self.cp_solver_name = cp_solver_name
        self.supported_cp_solvers = [
            'gecode', 'chuffed', 'cbc', 'gurobi', 'picat', 'scip', 'choco', 'or-tools']
        assert(self.cp_solver_name in self.supported_cp_solvers)
        self.cp_solver = minizinc.Solver.lookup(self.cp_solver_name)
        self.time_limit = time_limit
        self.cp_boolean_variables = []
        self.cp_constraints = ""
        self.mzn_file_name = f"warp_nr_{nrounds}.mzn"
        self.result_file_name = f"result_nr_{nrounds}.txt"
        self.permute_nibbles = [31, 6, 29, 14, 1, 12, 21, 8, 27, 2, 3, 0, 25, 4, 23, 10,
                                15, 22, 13, 30, 17, 28, 5, 24, 11, 18, 19, 16, 9, 20, 7, 26]
        # a0, b0: msb
        # convert monomial prediction table (MPT) to a boolean function
        # and then minimize its CNF (POS) representation via Quine-McCluskey algorithm
        self.sbox_predicate = "predicate sbox (var bool:a0, var bool:a1, var bool:a2, var bool:a3, var bool:b0, var bool:b1, var bool:b2, var bool:b3) = ((not a0) \\/ (not a3) \\/ b0 \\/ (not b1) \\/ (not b3)) /\\ ((not a0) \\/ a1 \\/ a3 \\/ b0 \\/ b1) /\\ (a1 \\/ a2 \\/ (not b2) \\/ (not b3)) /\\ (a2 \\/ (not b1) \\/ (not b3)) /\\ (a2 \\/ a3 \\/ (not b3)) /\\ ((not a0) \\/ (not a2) \\/ (not a3) \\/ (not b0) \\/ b1 \\/ (not b3)) /\\ (a0 \\/ a1 \\/ (not a3) \\/ b0 \\/ b1 \\/ b2) /\\ (a1 \\/ (not b1) \\/ (not b2)) /\\ (a1 \\/ a3 \\/ (not b2)) /\\ (a0 \\/ (not a2) \\/ a3 \\/ b3) /\\ ((not a2) \\/ (not b0) \\/ (not b1) \\/ b3) /\\ ((not a0) \\/ (not a1) \\/ (not a3) \\/ b0 \\/ b1 \\/ b3) /\\ ((not a1) \\/ (not b0) \\/ (not b1) \\/ b2) /\\ (a0 \\/ (not a1) \\/ a3 \\/ b2) /\\ ((not a2) \\/ b0 \\/ b1 \\/ b3) /\\ ((not a0) \\/ (not a1) \\/ (not a3) \\/ b2) /\\ ((not a0) \\/ (not a1) \\/ (not a2) \\/ (not b2) \\/ b3) /\\ ((not a1) \\/ a2 \\/ b0 \\/ b2 \\/ b3) /\\ ((not a1) \\/ (not a2) \\/ (not a3) \\/ b1 \\/ (not b2)) /\\ ((not a0) \\/ a1 \\/ a2 \\/ b1 \\/ b2 \\/ b3) /\\ ((not a1) \\/ (not a2) \\/ (not a3) \\/ b1 \\/ b3) /\\ (a2 \\/ (not a3) \\/ b1 \\/ b2 \\/ b3) /\\ (a0 \\/ a1 \\/ (not a2) \\/ (not b1) \\/ b3) /\\ ((not a1) \\/ a3 \\/ (not b0) \\/ b2 \\/ (not b3)) /\\ ((not a3) \\/ b0 \\/ (not b1) \\/ (not b2) \\/ (not b3)) /\\ (a1 \\/ (not b0) \\/ (not b2) \\/ (not b3)) /\\ ((not a1) \\/ a3 \\/ (not b1) \\/ b2 \\/ (not b3)) /\\ (a0 \\/ a1 \\/ a2 \\/ (not b3)) /\\ (a1 \\/ (not a2) \\/ a3 \\/ (not b1) \\/ b3);\n"

    def ordered_set(self, seq):
        """
        This method eliminates duplicated elements in a given list, 
        and returns a list in which each elements appears only once
        """

        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]

    @staticmethod
    def flatten_state(s):
            state_bits = [s[i][j] for i in range(len(s)) for j in range(len(s[0]))]
            return state_bits

    def inv_permute_nibbles(self, state):
        temp = [0]*32
        for i in range(32):
            temp[i] = state[self.permute_nibbles[i]]
        return temp    
    
    def generate_round_x_variables(self, rn):
        """
        Generate the input variables of rn'th round
        """

        x = [[f"x_{rn}[{nibble}, {bit}]" for bit in range(4)] for nibble in range(32)]
        self.cp_boolean_variables.extend(self.flatten_state(x))
        return x
    
    def generate_round_y_z_k_variables(self, rn):
        """
        Generate the intermediate variables in rn'th round
        """

        y = [[f"y_{rn}[{nibble}, {bit}]" for bit in range(4)] for nibble in range(16)]
        z = [[f"z_{rn}[{nibble}, {bit}]" for bit in range(4)] for nibble in range(16)]
        k = [[f"k[{rn}, {nibble}, {bit}]" for bit in range(4)] for nibble in range(16)]
        self.cp_boolean_variables.extend(self.flatten_state(y))
        self.cp_boolean_variables.extend(self.flatten_state(z))
        self.cp_boolean_variables.extend(self.flatten_state(k))
        return y, z, k

    @staticmethod
    def constraints_by_fork(a, b1, b2):
        """
        a ---fork---> (b1, b2)
        """
        constraints = f"constraint {a} <-> {b1} \\/ {b2};\n"
        return constraints
    
    def nibble_fork(self, s, s1, s2):
        constraints = ""
        for bit in range(4):
            constraints += self.constraints_by_fork(s[bit], s1[bit], s2[bit])
        return constraints

    @staticmethod
    def constraints_by_xor(a1, a2, b):
        """
        a1, a2 ---> b = a1 + a2
        """

        constraints = f"constraint {b} = {a1} + {a2};\n"
        return constraints
    
    def nibble_xor(self, s1, s2, s):
        constraints = ""
        for bit in range(4):
            constraints += self.constraints_by_xor(s1[bit], s2[bit], s[bit])
        return constraints
    
    @staticmethod
    def constraints_by_3xor(a1, a2, a3, b):
        """
        a1, a2, a3 ---> b = a1 + a2 + a3
        """

        constraints = f"constraint {b} = {a1} + {a2} + {a3};\n"
        return constraints
    
    def nibble_3xor(self, s1, s2, s3, s):
        constraints = ""
        for bit in range(4):
            constraints += self.constraints_by_3xor(s1[bit], s2[bit], s3[bit], s[bit])
        return constraints
    
    def constraints_by_sbox(self, variable1, variable2):
        """
        Generate the constraints by Sbox layer.
        """

        constraint = f"constraint sbox({variable1[0]}, {variable1[1]}, {variable1[2]}, {variable1[3]}, {variable2[0]}, {variable2[1]}, {variable2[2]}, {variable2[3]});\n"
        return constraint

    def generate_objective_function(self):
        """
        Create the objective function of the MILP model.
        """

        # objective = f"solve minimize {1};\n"
        objective = "solve satisfy;"
        return objective
    
    def generate_constraints(self):
        """
        Generate the constraints of MILP model
        """

        constraints = ""
        for rn in range(self.nrounds):
            x_in = self.generate_round_x_variables(rn)
            y, z, k = self.generate_round_y_z_k_variables(rn)
            x_out = self.generate_round_x_variables(rn + 1)
            x_middle = self.inv_permute_nibbles(x_out)
            for nibble in range(16):                
                constraints += self.nibble_fork(x_in[2*nibble], y[nibble], x_middle[2*nibble])
                constraints += self.constraints_by_sbox(y[nibble], z[nibble])
                constraints += self.nibble_3xor(z[nibble], k[nibble], x_in[2*nibble + 1], x_middle[2*nibble + 1])                
        return constraints
    
    def exclude_key_independent_term(self):
        """
        Limit the key variables to not be all zero
        """

        constraint = f"constraint exists(k) = 1;\n"
        return constraint
    
    def limit_the_output_choice_vector(self):
        """
        Limit the output choice vector to be unity
        """

        constraint = "constraint "
        x_out = self.generate_round_x_variables(self.nrounds)
        x_out_bits = self.flatten_state(x_out)
        constraint += " + ".join(x_out_bits) + " <= 1;\n"
        constraint = f"sum(x[{self.nrounds}, 0..31, 0..3]) <= 1;\n"
        return constraint
        
    def declare_boolean_vars(self):
        """
        Declare binary variables of MILP model
        """

        boolean_variables = ""
        for rn in range(self.nrounds + 1):
            boolean_variables += f"array [0..31, 0..3] of var bool: x_{rn};\n"
        for rn in range(self.nrounds):
            boolean_variables += f"array [0..15, 0..3] of var bool: y_{rn};\n"
            boolean_variables += f"array [0..15, 0..3] of var bool: z_{rn};\n"
        boolean_variables += f"array [0..{self.nrounds - 1}, 0..15, 0..3] of var bool: k;\n"
        return boolean_variables
    

    def generate_mzn_contents(self):
        """
        Generate the MILP model describing the propagation monomial prediction vectors
        """

        mzn_contents = "% Integral attack on {} rounds of WARP\n".format(self.nrounds)
        mzn_contents += self.declare_boolean_vars()
        mzn_contents += self.sbox_predicate
        mzn_contents += self.generate_constraints()
        mzn_contents += self.exclude_key_independent_term()
        return mzn_contents

    def check_cube(self, fixed_indices=[0]):
        mzn_file_contents_main_body = self.generate_mzn_contents()
        # Fix the input choice vector
        input_active_pattern = [str(not(i in fixed_indices)).lower() for i in range(128)]
        input_active_pattern = ", ".join(input_active_pattern)
        input_active_pattern = f"x_0 = array2d(0..31, 0..3, [{input_active_pattern}]);\n"
        
        if self.time_limit != -1:
            time_limit = datetime.timedelta(seconds=self.time_limit)
        else:
            time_limit = None

        balanced_bits = []
        not_checked_bits = []
        start_time = time.time()
        for output_bit in range(128):
            output_active_pattern = [str(i == output_bit).lower() for i in range(128)]
            output_active_pattern = ", ".join(output_active_pattern)
            output_active_pattern = f"x_{self.nrounds} = array2d(0..31, 0..3, [{output_active_pattern}]);\n"
            boundary_conditions = input_active_pattern + output_active_pattern
            mzn_file_contents = mzn_file_contents_main_body + boundary_conditions + self.generate_objective_function()
            with open(self.mzn_file_name, "w") as mzn_file:
                mzn_file.write(mzn_file_contents)
            ##########################
            ##########################
            self.cp_model = minizinc.Model()
            self.cp_model.add_file(self.mzn_file_name)
            self.cp_inst = minizinc.Instance(solver=self.cp_solver, model=self.cp_model)
            result = self.cp_inst.solve(timeout=time_limit)
            ##########################
            ##########################
            if result.status == minizinc.Status.OPTIMAL_SOLUTION or result.status == minizinc.Status.SATISFIED or \
                                result.status == minizinc.Status.ALL_SOLUTIONS:       
                print("Output bit number {:03d} may NOT be key-independent :-(".format(output_bit))
            elif result.status == minizinc.Status.UNSATISFIABLE:
                balanced_bits.append(output_bit)
                print("Output bit number {:03d} is key-independent ;-)".format(output_bit))
            else:
                not_checked_bits.append(output_bit)
                print("Output bit number {:03d} was not checked due to the time limit".format(output_bit))
        os.remove(self.mzn_file_name)
        elapsed_time = time.time() - start_time
        number_of_balanced_bits = len(balanced_bits)
        print(f"Number of key-independent bits: {number_of_balanced_bits}")
        print(f"Key-Independent bits:\n{balanced_bits}")
        print(f"Not-Checked bits:{not_checked_bits}\n")
        print("Time used to solve: {:0.02f}".format(elapsed_time))        
        ######################### Save results in output file ##############################
        with open(self.result_file_name, "a") as outputfile:
            separator = "#"*100 + "\n"
            outputfile.write(separator)
            outputfile.write(f"Fixed input positions: {fixed_indices}\n")
            outputfile.write(f"Key-independent output positions: {balanced_bits}\n")
            outputfile.write(f"Number of key-independent bits: {number_of_balanced_bits}\n")
        ####################################################################################
        return number_of_balanced_bits

def parse_args():
    """
    parse input parameters
    """

    parser = ArgumentParser(description="This tool derives and solves the CP "
                                        "model corresponding to integral analysis "
                                        "based on monomial prediction",
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-nr", "--nrounds", default=21, type=int, help="number of rounds\n")
    parser.add_argument("-sl", "--solver", default="chuffed", type=str,
                        choices=['gecode', 'chuffed', 'coin-bc', 'gurobi', 'picat', 'scip', 'choco', 'or-tools'],
                        help="choose a cp solver\n")
    parser.add_argument("-tl", "--timelimit", default=5, type=int, help="set a time limit for the solver in seconds\n")
    return vars(parser.parse_args())
    
if __name__ == '__main__':
    locals().update(parse_args())
    print(nrounds, solver, timelimit)
    warp = Warp(nrounds=nrounds,
                cp_solver_name=solver,
                time_limit=timelimit)
    with open(warp.result_file_name, "w") as outputfile:
        outputfile.write(f"Results of applying monomial prediction method on {warp.nrounds} rounds of WARP\n")
    nb_dict = dict()
    warp.generate_mzn_contents()
    for bit in range(0, 128):
        separator = "#"*100 + "\n"
        print(separator)
        print(f"Fixed input index: {bit}")
        number_of_balanced_bits = warp.check_cube(fixed_indices=[bit])
        nb_dict[bit] = number_of_balanced_bits
    print(nb_dict)
