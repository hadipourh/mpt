"""
Applying the monomial prediction technique to find integral
distinguishers of WARP in the single-key setting
Copyright (C) Jan 2, 2022 Hosein Hadipour

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
import pickle
from itertools import combinations
from math import factorial
from math import log2
from typing import final
from pysat import solvers
from pysat import formula
from argparse import ArgumentParser, RawTextHelpFormatter
from rich.progress import Progress, TimeElapsedColumn, BarColumn
# from pysat import pb
# import datetime

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
    def __init__(self, nrounds=10, solver_name=solvers.SolverNames.cadical153):
        Warp.count += 1
        self.nrounds = nrounds
        self.sat_solver_name = solver_name
        self.supported_sat_solvers = [solver for solver in solvers.SolverNames.__dict__.keys() if not solver.startswith('__')]
        if self.sat_solver_name not in self.supported_sat_solvers:
            raise ValueError(f"Unsupported SAT solver: {self.sat_solver_name}")
        self.sat_solver = solvers.Solver(name=self.sat_solver_name)
        self.cnf_formula = formula.CNF()
        self.solver = solvers.Solver(name=self.sat_solver_name)
        self.variables_dictionary = dict()
        self.top_variable_identifier_so_far = 0
        self.cnf_file_name = f"warp_nr_{self.nrounds}.cnf"
        self.result_file_name = f"result_nr_{self.nrounds}.txt"
        self.permute_nibbles = [31, 6, 29, 14, 1, 12, 21, 8, 27, 2, 3, 0, 25, 4, 23, 10,
                                15, 22, 13, 30, 17, 28, 5, 24, 11, 18, 19, 16, 9, 20, 7, 26]
        self.RC0 = [0x0, 0x0, 0x1, 0x3, 0x7, 0xf, 0xf, 0xf, 0xe, 0xd, 0xa, 0x5, 0xa, 0x5, 0xb, 0x6, 0xc, 0x9, 0x3, 0x6, 0xd, 0xb, 0x7, 0xe, 0xd, 0xb, 0x6, 0xd, 0xa, 0x4, 0x9, 0x2, 0x4, 0x9, 0x3, 0x7, 0xe, 0xc, 0x8, 0x1, 0x2]
        self.RC1 = [0x4, 0xc, 0xc, 0xc, 0xc, 0xc, 0x8, 0x4, 0x8, 0x4, 0x8, 0x4, 0xc, 0x8, 0x0, 0x4, 0xc, 0x8, 0x4, 0xc, 0xc, 0x8, 0x4, 0xc, 0x8, 0x4, 0x8, 0x0, 0x4, 0x8, 0x0, 0x4, 0xc, 0xc, 0x8, 0x0, 0x0, 0x4, 0x8, 0x4, 0xc]
        # 1: msb of input
        # 5: msb of output
        # input bits: 1, 2, 3, 4
        # output bits: 5, 6, 7, 8
        # convert monomial prediction table (MPT) to a boolean function
        # and then minimize its CNF (POS) representation via Quine-McCluskey algorithm
        self.sbox_cnf_pattern = [[-1, -4, 5, -6, -8],
                                [-1, 2, 4, 5, 6],
                                [2, 3, -7, -8],
                                [3, -6, -8],
                                [3, 4, -8],
                                [-1, -3, -4, -5, 6, -8],
                                [1, 2, -4, 5, 6, 7],
                                [2, -6, -7],
                                [2, 4, -7],
                                [1, -3, 4, 8],
                                [-3, -5, -6, 8],
                                [-1, -2, -4, 5, 6, 8],
                                [-2, -5, -6, 7],
                                [1, -2, 4, 7],
                                [-3, 5, 6, 8],
                                [-1, -2, -4, 7],
                                [-1, -2, -3, -7, 8],
                                [-2, 3, 5, 7, 8],
                                [-2, -3, -4, 6, -7],
                                [-1, 2, 3, 6, 7, 8],
                                [-2, -3, -4, 6, 8],
                                [3, -4, 6, 7, 8],
                                [1, 2, -3, -6, 8],
                                [-2, 4, -5, 7, -8],
                                [-4, 5, -6, -7, -8],
                                [2, -5, -7, -8],
                                [-2, 4, -6, 7, -8],
                                [1, 2, 3, -8],
                                [2, -3, 4, -6, 8]]
        # xor inputs: 1, 2
        # xor output: 3
        self.xor_cnf_template = [[-1, -2], [-2, 3], [-1, 3], [1, 2, -3]]

        # xor3 inputs: 1, 2, 3
        # xor3 output: 4
        self.xor3_cnf_template = [[-2, -3], [-1, -3], [-1, -2], [-3, 4], [-2, 4], [-1, 4], [1, 2, 3, -4]]
        self.xor3_with_constant_cnf_template = [[-2, -3], [-1, -3], [-1, -2], [-3, 4], [-2, 4], [-1, 4]]
    
    def generate_sbox_constraint(self, input_bits, output_bits):
        substitution_list = [self.variables_dictionary[x] for x in input_bits + output_bits]
        for sl in self.sbox_cnf_pattern:
            temp = []
            for index in sl:
                if index > 0:
                    temp.append(substitution_list[index - 1])
                else:
                    temp.append(-substitution_list[abs(index) - 1])
            self.cnf_formula.append(temp)

    def update_variables_dictionary(self, new_vars):
        """
        This method is used to update variables' dictionary
        """
        for nv in new_vars:
            if nv not in self.variables_dictionary.keys():
                self.top_variable_identifier_so_far += 1
                self.variables_dictionary[nv] = self.top_variable_identifier_so_far

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

        x = [[f"x_{rn}_{nibble}_{bit}" for bit in range(4)] for nibble in range(32)]
        self.update_variables_dictionary(self.flatten_state(x))
        return x
    
    def generate_round_y_z_k_variables(self, rn):
        """
        Generate the intermediate variables in rn'th round
        """

        y = [[f"y_{rn}_{nibble}_{bit}" for bit in range(4)] for nibble in range(16)]
        z = [[f"z_{rn}_{nibble}_{bit}" for bit in range(4)] for nibble in range(16)]
        k = [[f"k_{rn}_{nibble}, {bit}" for bit in range(4)] for nibble in range(16)]
        self.update_variables_dictionary(self.flatten_state(y))
        self.update_variables_dictionary(self.flatten_state(z))
        self.update_variables_dictionary(self.flatten_state(k))
        return y, z, k
        
    def constraints_by_fork(self, a, b1, b2):
        """
        a ---fork---> (b1, b2)
        """

        self.cnf_formula.append([-self.variables_dictionary[a],\
                                self.variables_dictionary[b1],\
                                self.variables_dictionary[b2]])
        self.cnf_formula.append([-self.variables_dictionary[b1], self.variables_dictionary[a]])
        self.cnf_formula.append([-self.variables_dictionary[b2], self.variables_dictionary[a]])
    
    def nibble_fork(self, s, s1, s2):
        constraints = ""
        for bit in range(4):
            self.constraints_by_fork(s[bit], s1[bit], s2[bit])
        return constraints

    def constraints_by_xor(self, a1, a2, b):
        """
        a1, a2 ---> b = a1 + a2
        """
        substitution_list = [self.variables_dictionary[x] for x in  [a1, a2, b]]
        for sl in self.xor_cnf_template:
            temp = []
            for index in sl:
                if index > 0:
                    temp.append(substitution_list[index - 1])
                else:
                    temp.append(-substitution_list[abs(index) - 1])
            self.cnf_formula.append(temp)
    
    def nibble_xor(self, s1, s2, s):
        for bit in range(4):
            self.constraints_by_xor(s1[bit], s2[bit], s[bit])
    
    def constraints_by_3xor(self, a1, a2, a3, b, constant=0):
        """
        a1, a2, a3 ---> b = a1 + a2 + a3
        or
        a1, a2, a3 ---> b = a1 + a2 + a3 + 1
        """
        
        assert(constant in [0, 1])
        substitution_list = [self.variables_dictionary[x] for x in  [a1, a2, a3, b]]        
        if constant == 0:
            template = self.xor3_cnf_template
        elif constant == 1:
            template = self.xor3_with_constant_cnf_template        
        for sl in template:
            temp = []
            for index in sl:
                if index > 0:
                    temp.append(substitution_list[index - 1])
                else:
                    temp.append(-substitution_list[abs(index) - 1])
            self.cnf_formula.append(temp)
    
    def nibble_3xor(self, s1, s2, s3, s, constant):
        for bit in range(4):
            constant_bit = (constant >> (3 - bit)) & 0x1
            self.constraints_by_3xor(s1[bit], s2[bit], s3[bit], s[bit], constant_bit)
    
    def generate_constraints(self):
        """
        Generate the constraints of MILP model
        """

        for rn in range(self.nrounds):
            x_in = self.generate_round_x_variables(rn)
            y, z, k = self.generate_round_y_z_k_variables(rn)
            x_out = self.generate_round_x_variables(rn + 1)
            x_middle = self.inv_permute_nibbles(x_out)
            for nibble in range(16):                
                self.nibble_fork(x_in[2*nibble], y[nibble], x_middle[2*nibble])
                self.generate_sbox_constraint(y[nibble], z[nibble])                    
                if nibble == 1:
                    constant = self.RC0[rn]
                elif nibble == 3:
                    constant = self.RC1[rn]
                else:
                    constant = 0
                self.nibble_3xor(z[nibble], k[nibble], x_in[2*nibble + 1], x_middle[2*nibble + 1], constant)
    
    def exclude_key_independent_term(self):
        """
        Limit the key variables to not be all zero
        """

        key_vars = [self.variables_dictionary[var] \
            for var in self.variables_dictionary.keys() if "k" in var]
        self.cnf_formula.append(key_vars)
    
    def generate_sat_model(self):
        self.generate_constraints()
        self.exclude_key_independent_term()
        self.cnf_formula.to_file(self.cnf_file_name)
        self.sat_solver.append_formula(self.cnf_formula)

    def check_cube(self, fixed_indices=[0], target_output_bits = range(128)):
        # Fix the input choice vector
        input_vars = self.flatten_state(self.generate_round_x_variables(0))
        output_vars = self.flatten_state(self.generate_round_x_variables(self.nrounds))

        input_active_pattern = []
        for i in range(128):
            if i in fixed_indices:
                input_active_pattern.append(-self.variables_dictionary[input_vars[i]])
            else:
                input_active_pattern.append(self.variables_dictionary[input_vars[i]])
        
        balanced_bits = []
        not_checked_bits = []
        start_time = time.time()
        with Progress("[progress.description]{task.description}", 
                  BarColumn(),  # Add a bar to visually show progress
                  "[progress.percentage]{task.percentage:>3.0f}%", 
                  "•", 
                  TimeElapsedColumn()) as progress:
            task = progress.add_task("Checking output bits...", total=len(target_output_bits))        
            for output_bit in target_output_bits:
                output_active_pattern = []
                for i in range(128):
                    if i != output_bit:
                        output_active_pattern.append(-self.variables_dictionary[output_vars[i]])
                    else:
                        output_active_pattern.append(self.variables_dictionary[output_vars[i]])
                assumptions = input_active_pattern + output_active_pattern
                ##########################
                ##########################
                result = self.sat_solver.solve(assumptions=assumptions)
                ##########################
                ##########################
                if result == True:
                    # print("Output bit number {:03d} may NOT be key-independent :-(".format(output_bit))
                    pass
                elif result == False:
                    balanced_bits.append(output_bit)
                    # print("Output bit number {:03d} is key-independent ;-)".format(output_bit))
                else:
                    not_checked_bits.append(output_bit)
                    # print("Output bit number {:03d} was not checked!".format(output_bit))
                progress.update(task, advance=1, elapsed_time=f"{time.time() - start_time:0.02f}s")
        elapsed_time = time.time() - start_time
        number_of_balanced_bits = len(balanced_bits)
        print(f"Number of key-independent bits: {number_of_balanced_bits}")
        print(f"Key-Independent bits:\n{balanced_bits}")
        print(f"Not-Checked bits:{not_checked_bits}\n")
        print("Time used to solve: {:0.02f}".format(elapsed_time))        
        ######################### Save results in output file ##############################
        with open(self.result_file_name, "a") as outputfile:
            separator_line = "#"*100 + "\n"
            outputfile.write(separator_line)
            outputfile.write(f"Fixed input positions: {fixed_indices}\n")
            outputfile.write(f"Key-independent output positions: {balanced_bits}\n")
            outputfile.write(f"Number of key-independent bits: {number_of_balanced_bits}\n")
        ####################################################################################
        return balanced_bits

def parse_args():
    """
    parse input parameters
    """

    parser = ArgumentParser(description="This tool derives and solves the SAT "
                                        "model corresponding to integral analysis "
                                        "based on monomial prediction",
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-nr", "--nrounds", default=21, type=int, help="number of rounds\n")
    parser.add_argument("-sl", "--solver", default="minisat22", type=str,
                        choices=['cadical', 'glucose3', 'glucose4', 'lingeling', 'maplechrono', 'maplecm', 'maplesat', 'minicard', 'minisat22', 'minisat-gh'],
                        help="choose a SAT solver\n")
    return vars(parser.parse_args())

if __name__ == '__main__':
    locals().update(parse_args())    
    separator_line = "#"*100 + "\n"
    warp = Warp(nrounds=nrounds,
                solver_name=solver)
    with open(warp.result_file_name, "w") as outputfile:
        outputfile.write(f"Results of applying monomial prediction method to {warp.nrounds} rounds of WARP\n")
    mp_dict = dict()
    warp.generate_sat_model()
    if not os.path.exists("mpdict.pyobj"):
        start_time = time.time()
        for bit in range(0, 128):            
            print(separator_line)
            print(f"Fixed index: {bit}")
            mp_dict[(bit,)] = warp.check_cube(fixed_indices=[bit])
        elapsed_time = time.time() - start_time
        time_line = "Total time to compute monomial prediction dictionary (mpdict): {:0.02f}\n".format(elapsed_time)        
        with open(warp.result_file_name, "a") as outputfile:
            outputfile.write(time_line)
        print(time_line)
        # Store the results in non-volatile memory for later uses
        with open("mpdict.pyobj", "wb") as mp_dict_file:
            pickle.dump(mp_dict, mp_dict_file)
    else:
        with open("mpdict.pyobj", "rb") as mp_dict_file:
            mp_dict = pickle.load(mp_dict_file)

    # Reduction of data complexity
    start_time = time.time()
    sufficient_list = [i for i in range(128) if mp_dict[(i,)] != []]
    sf_size = len(sufficient_list)
    constant_part_size = 2
    flag = True
    while flag:
        flag = False
        print(separator_line)
        print(f"Checking constant parts of size {constant_part_size} bits")
        total_possible_cases = factorial(sf_size) // \
                         (factorial(sf_size - constant_part_size) * factorial(constant_part_size))
        total_possible_cases = log2(total_possible_cases)
        print("Worst-case complexity: 2^({:0.02f})".format(total_possible_cases))        
        time.sleep(5)
        for fixed_indices in combinations(sufficient_list, constant_part_size):
            balanced_sets = [set(mp_dict[(i,)]) for i in fixed_indices]
            common_output_bbits = set.intersection(*balanced_sets)
            if common_output_bbits != {}:
                key_independent_bits = warp.check_cube(fixed_indices=list(fixed_indices), 
                                                        target_output_bits=list(common_output_bbits))
                if key_independent_bits != []:
                    flag = True
                    constant_part_size += 1
                    break
    
    # Find a set of fixed inidces yielding the maximum number of key-independent bits
    best_fixed_positions = (128,)
    mp_dict[(128,)] = []
    constant_part_size -= 1
    key_independent_bits = []
    with open(warp.result_file_name, "a") as outputfile:
        for fixed_indices in combinations(sufficient_list, constant_part_size):
            balanced_sets = [set(mp_dict[(i,)]) for i in fixed_indices]
            common_output_bbits = set.intersection(*balanced_sets)
            if common_output_bbits != {}:
                key_independent_bits = warp.check_cube(fixed_indices=list(fixed_indices), 
                                                            target_output_bits=list(common_output_bbits))
                if len(key_independent_bits) >= 1:
                    outputfile.write(separator_line)
                    outputfile.write(f"{fixed_indices} ==> {key_independent_bits}\n")
                if len(key_independent_bits) > len(mp_dict[best_fixed_positions]):
                    best_fixed_positions = fixed_indices
                    mp_dict[best_fixed_positions] = key_independent_bits
    
    elapsed_time = time.time() - start_time
    time_line = "Time used to reduce the data complexity: {:0.02f}\n".format(elapsed_time)
    print(time_line)
    final_result = f"Best fixed input positions: {best_fixed_positions}\n"
    final_result += f"Key-independent output positions: {mp_dict[tuple(best_fixed_positions)]}"
    print(final_result)
    with open(warp.result_file_name, "a") as output_file:
        output_file.write(separator_line)
        output_file.write(time_line)
        output_file.write(final_result)
