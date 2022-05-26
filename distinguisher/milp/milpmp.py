"""
Applying the monomial prediction technique to find integral
distinguishers of WARP in the single-key setting
Copyright (C) Dec, 28, 2021  Hosein Hadipour

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

from email.policy import default
from gurobipy import *
import os
import time
from argparse import ArgumentParser, RawTextHelpFormatter
from matplotlib.cbook import flatten

from plumbum import local


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
    def __init__(self, nrounds=10):
        Warp.count += 1
        self.nrounds = nrounds
        self.milp_variables = []
        self.lp_file_name = f"warp_nr_{nrounds}.lp"
        self.result_file_name = f"result_nr_{nrounds}.txt"
        self.permute_nibbles = [31, 6, 29, 14, 1, 12, 21, 8, 27, 2, 3, 0, 25, 4, 23, 10,
                                15, 22, 13, 30, 17, 28, 5, 24, 11, 18, 19, 16, 9, 20, 7, 26]
        # a0, b0: msb
        # Method 1: convert monomial prediction table (MPT) to a boolean function
        # and then minimize its CNF (POS) represntation via Quine–McCluskey algorithm
        # ***This mehod yields a much better performance when we use method 1 to check cubes***
        self.sbox_mpt_ineqs = ["- a0 - a3 + b0 - b1 - b3 >= -3",
                                "- a0 + a1 + a3 + b0 + b1 >= 0",
                                "a1 + a2 - b2 - b3 >= -1",
                                "a2 - b1 - b3 >= -1",
                                "a2 + a3 - b3 >= 0",
                                "- a0 - a2 - a3 - b0 + b1 - b3 >= -4",
                                "a0 + a1 - a3 + b0 + b1 + b2 >= 0",
                                "a1 - b1 - b2 >= -1",
                                "a1 + a3 - b2 >= 0",
                                "a0 - a2 + a3 + b3 >= 0",
                                "- a2 - b0 - b1 + b3 >= -2",
                                "- a0 - a1 - a3 + b0 + b1 + b3 >= -2",
                                "- a1 - b0 - b1 + b2 >= -2",
                                "a0 - a1 + a3 + b2 >= 0",
                                "- a2 + b0 + b1 + b3 >= 0",
                                "- a0 - a1 - a3 + b2 >= -2",
                                "- a0 - a1 - a2 - b2 + b3 >= -3",
                                "- a1 + a2 + b0 + b2 + b3 >= 0",
                                "- a1 - a2 - a3 + b1 - b2 >= -3",
                                "- a0 + a1 + a2 + b1 + b2 + b3 >= 0",
                                "- a1 - a2 - a3 + b1 + b3 >= -2",
                                "a2 - a3 + b1 + b2 + b3 >= 0",
                                "a0 + a1 - a2 - b1 + b3 >= -1",
                                "- a1 + a3 - b0 + b2 - b3 >= -2",
                                "- a3 + b0 - b1 - b2 - b3 >= -3",
                                "a1 - b0 - b2 - b3 >= -2",
                                "- a1 + a3 - b1 + b2 - b3 >= -2",
                                "a0 + a1 + a2 - b3 >= 0",
                                "a1 - a2 + a3 - b1 + b3 >= -1"]

        # Method 2: Derive the convex-hull of valid monomial trails and 
        # then reduce the number of inequalities according to a heuristic algorithm
        # self.sbox_mpt_ineqs = ["-2 a0 - a1 - a2 - a3 - b0 + b1 + b2 - b3 >= -5",
        #                         "-2 a0 + a1 + a2 + 2 a3 + b0 + 2 b1 + b2 >= 0",
        #                         "-a0 - 2 a1 - b0 - b1 + 2 b2 - b3 >= -4",
        #                         "-a0 - a1 - a3 + 2 b0 + 2 b1 + 2 b2 + 2 b3 >= 0",
        #                         "-a0 - a1 + a2 - b1 + b2 - b3 >= -2",
        #                         "-a0 - a2 - 2 a3 + 2 b0 - 2 b1 - b2 - b3 >= -6",
        #                         "-2 a1 - 2 a2 - a3 + 2 b1 - b2 + b3 >= -4",
        #                         "-a1 + a2 - 2 a3 + b1 + 2 b2 + 2 b3 >= -1",
        #                         "-a2 - b0 - b1 + b3 >= -2",
        #                         "a0 + a1 + a2 - b3 >= 0",
        #                         "a1 - a2 + a3 + b0 - b2 + b3 >= 0",
        #                         "a1 - b0 - b2 - b3 >= -2",
        #                         "a1 + a2 - b1 - b2 - b3 >= -1",
        #                         "a1 + a2 + a3 - b2 - b3 >= 0",
        #                         "a0 - a1 + a3 + b2 >= 0",
        #                         "a0 + a1 - 2 a2 + a3 - b1 + 2 b3 >= -1",
        #                         "a0 + a1 - a3 + b0 + b1 + b2 >= 0"]

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

        x = [[f"x_{rn}_{nibble}_{bit}" for bit in range(4)] for nibble in range(32)]
        self.milp_variables.extend(self.flatten_state(x))
        return x
    
    def generate_round_y_z_k_variables(self, rn):
        """
        Generate the intermediate variables in rn'th round
        """

        y = [[f"y_{rn}_{nibble}_{bit}" for bit in range(4)] for nibble in range(16)]
        z = [[f"z_{rn}_{nibble}_{bit}" for bit in range(4)] for nibble in range(16)]
        k = [[f"sk_{rn}_{nibble}_{bit}" for bit in range(4)] for nibble in range(16)]
        self.milp_variables.extend(self.flatten_state(y))
        self.milp_variables.extend(self.flatten_state(z))
        self.milp_variables.extend(self.flatten_state(k))
        return y, z, k

    @staticmethod
    def constraints_by_fork(a, b1, b2):
        """
        a ---fork---> (b1, b2)
        """

        constraints = f"{b1} + {b2} - {a} >= 0\n"
        constraints += f"{a} - {b1} >= 0\n"
        constraints += f"{a} - {b2} >= 0\n" 
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

        constraints = f"{a1} + {a2} - {b} = 0\n"
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

        constraints = f"{a1} + {a2} + {a3} - {b} = 0\n"
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

        constraints = ""
        for ineq in self.sbox_mpt_ineqs:
            temp = ineq
            for i in range(4):
                temp = temp.replace("a{}".format(i), variable1[i])
            for i in range(4):
                temp = temp.replace("b{}".format(i), variable2[i])
            constraints += temp + "\n"
        return constraints
    
    def generate_key_schedule_constraints(self):
        """
        Modelize the key schedule of WARP
        """

        constraints = "\ngeneral constraints\n"
        k0 = [f"k0_{bit}" for bit in range(64)]
        k1 = [f"k1_{bit}" for bit in range(64)]
        self.milp_variables.extend(k0 + k1)
        for bit_position in range(64):
            nibble = bit_position // 4
            bit = bit_position % 4
            sub_keys_k0 = [f"sk_{rn}_{nibble}_{bit}" for rn in range(self.nrounds) if rn % 2 == 0]
            sub_keys_k0 = " , ".join(sub_keys_k0)            
            constraints += f"{k0[bit_position]} = OR ( {sub_keys_k0} )\n"
            sub_keys_k1 = [f"sk_{rn}_{nibble}_{bit}" for rn in range(self.nrounds) if rn % 2 == 1]
            sub_keys_k1 = " , ".join(sub_keys_k1)
            constraints += f"{k1[bit_position]} = OR ( {sub_keys_k1} )\n"
        return constraints

    def generate_objective_function(self):
        """
        Create the objective function of the MILP model.
        """

        objective_function = "Minimize\n"        
        objective_function += "1\n"
        return objective_function
    
    def generate_constraints(self):
        """
        Generate the constraints of MILP model
        """
        
        constraints = "subject to\n"
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

        key_vars = [v for v in self.milp_variables if "k" in v]
        constraint = " + ".join(key_vars) + " >= 100\n"
        return constraint
    
    def limit_the_output_choice_vector(self):
        """
        Limit the output choice vector to be unity
        """

        x_out = self.generate_round_x_variables(self.nrounds)
        x_out_bits = self.flatten_state(x_out)
        constraint = " + ".join(x_out_bits) + " <= 1\n"
        return constraint
        
    def declare_binary_vars(self):
        """
        Declare binary variables of MILP model
        """
        
        self.milp_variables = self.ordered_set(self.milp_variables)
        constraints = "Binary\n"
        constraints += "\n".join(self.milp_variables)
        return constraints
    

    def make_model(self):
        """
        Generate the MILP model describing the propagation monomial prediction vectors
        """

        lp_contents = "\\ Integral attack on {} rounds of WARP\n".format(self.nrounds)
        lp_contents += self.generate_objective_function()
        lp_contents += self.generate_constraints()
        lp_contents += self.exclude_key_independent_term()     
        # lp_contents += self.limit_the_output_choice_vector()
        lp_contents += self.declare_binary_vars()
        lp_contents += self.generate_key_schedule_constraints()
        lp_contents += "End"
        with open(self.lp_file_name, "w") as lp_file:
            lp_file.write(lp_contents)
        
    def check_cube_method1(self, fixed_indices=[0]):
        self.make_model()
        milp_model = read(self.lp_file_name)
        # os.remove(self.lp_file_name)
        milp_model.setParam(GRB.Param.OutputFlag, False)
        # milp_model.setParam(GRB.Param.Presolve, 1)
        input_mask_vars = self.flatten_state(self.generate_round_x_variables(0))
        input_mask_vars = [milp_model.getVarByName(x) for x in input_mask_vars]
        output_mask_vars = self.flatten_state(self.generate_round_x_variables(self.nrounds))
        output_mask_vars = [milp_model.getVarByName(x) for x in output_mask_vars]

        # Fix the input choice vector
        for i in range(128):
            if i in fixed_indices:
                milp_model.addConstr(input_mask_vars[i] == 0)
            else:
                milp_model.addConstr(input_mask_vars[i] == 1)

        # Limit the hamming weight of the output choice vector to be at most 1
        output_constraint = milp_model.addConstr((sum(output_mask_vars[i] for i in range(128)) <= 1))

        balanced_bits = list(range(128))
        start_time = time.time()
        milp_model.optimize()
        counter = 0
        while(milp_model.Status == GRB.Status.OPTIMAL):            
            for output_bit in range(128):
                if int(output_mask_vars[output_bit].Xn) == 1:
                    balanced_bits.remove(output_bit)
                    milp_model.addConstr(output_mask_vars[output_bit] == 0)
                    print("{:03d}: Output bit number {:03d} may NOT be key-independent :-(".format(counter, output_bit))
                    break
            counter += 1
            milp_model.optimize()
        elapsed_time = time.time() - start_time
        number_of_balanced_bits = len(balanced_bits)
        print(f"Number of balanced bits: {number_of_balanced_bits}")
        print(f"Balanced bits:\n{balanced_bits}")
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

    def check_cube_method2(self, fixed_indices=[0], time_limit=5):
        self.make_model()
        milp_model = read(self.lp_file_name)
        # os.remove(self.lp_file_name)
        milp_model.setParam(GRB.Param.OutputFlag, False)
        # milp_model.setParam(GRB.Param.Presolve, 0)
        milp_model.Params.TIME_LIMIT = time_limit
        input_mask_vars = self.flatten_state(self.generate_round_x_variables(0))
        input_mask_vars = [milp_model.getVarByName(x) for x in input_mask_vars]
        output_mask_vars = self.flatten_state(self.generate_round_x_variables(self.nrounds))
        output_mask_vars = [milp_model.getVarByName(x) for x in output_mask_vars]
        # Fix the input choice vector
        for i in range(128):
            if i in fixed_indices:
                milp_model.addConstr(input_mask_vars[i] == 0)
            else:
                milp_model.addConstr(input_mask_vars[i] == 1)

        balanced_bits = []
        not_checked_bits = []
        start_time = time.time()
        for output_bit in range(128):
            temporary_constraint1 = milp_model.addConstrs(
                (output_mask_vars[bit] == 0 for bit in range(128) if bit != output_bit), 
                name='output_constraint1')
            temporary_constraint2 = milp_model.addConstr(output_mask_vars[output_bit] == 1,
                name='output_constraint2')
            milp_model.optimize()
            if (milp_model.Status == GRB.Status.INFEASIBLE):
                balanced_bits.append(output_bit)
                print("Output bit number {:03d} is key-independent ;-)".format(output_bit))
            elif (milp_model.Status == GRB.Status.OPTIMAL):
                print("Output bit number {:03d} may NOT be key-independent :-(".format(output_bit))
            elif (milp_model.Status == GRB.TIME_LIMIT or milp_model.Status == GRB.INTERRUPTED):
                not_checked_bits.append(output_bit)
                print("Output bit number {:03d} was not checked due to the time limit".format(output_bit))
            milp_model.remove(temporary_constraint1)
            milp_model.remove(temporary_constraint2)
            milp_model.update()            
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

    def superpoly_recovery(self, active_indices=[24, 25, 26, 27], target_output_bit=22):
        """
        WARNING!
        The Gurobi PoolSearchMode supports at most 2,000,000,000
        solutions. As the superpoly becomes increasingly complicated, the actual number
        of trails is likely to surpass this bound.

        In this method, we set the cube variables to 1 and leave the non-cube variables
        as well as key variables to be free. Next, we find all possible solutions satisfying the
        propagation rules for monomial prediction. Lastly, we partition the solutions according to
        the key part and determine the parity of each subset. The key monomials whose corresponding subset of 
        solutions have odd entries are present in ANF of the targeted output bit.
        """

        fixed_indices = [i for i in range(128) if i not in active_indices]
        self.make_model()
        milp_model = read(self.lp_file_name)
        milp_model.Params.PoolSearchMode = 2
        milp_model.Params.PoolSolutions = GRB.MAXINT # MAXINT = 2000000000

        input_mask_vars = self.flatten_state(self.generate_round_x_variables(0))
        input_mask_vars = [milp_model.getVarByName(x) for x in input_mask_vars]
        output_mask_vars = self.flatten_state(self.generate_round_x_variables(self.nrounds))
        output_mask_vars = [milp_model.getVarByName(x) for x in output_mask_vars]
        non_cube_variables = [input_mask_vars[i] for i in fixed_indices]
        key_variables_0 = [f"k0_{bit}" for bit in range(64)]
        key_variables_0 = [milp_model.getVarByName(x) for x in key_variables_0]
        key_variables_1 = [f"k1_{bit}" for bit in range(64)]
        key_variables_1 = [milp_model.getVarByName(x) for x in key_variables_1]
        target_variables = non_cube_variables + key_variables_0 + key_variables_1
        # Set the cube variables to 1
        for i in range(128):
            if i not in fixed_indices:
                milp_model.addConstr(input_mask_vars[i] == 1)
        
        # Fixed output choice vector
        temporary_constraint1 = milp_model.addConstrs(
            (output_mask_vars[bit] == 0 for bit in range(128) if bit != target_output_bit), 
            name='output_constraint1')
        temporary_constraint2 = milp_model.addConstr(output_mask_vars[target_output_bit] == 1,
            name='output_constraint2')
        ###############################
        ###############################
        start_time = time.time()
        milp_model.optimize()
        number_of_solutions = milp_model.SolCount
        elapsed_time = time.time() - start_time
        ###############################
        ###############################
        print("Time used to solve the MILP problem: {:.02f}".format(elapsed_time))
        print("Number of solutions: {}".format(number_of_solutions)) 

        support = dict()
        for sn in range(number_of_solutions):
            milp_model.Params.SolutionNumber = sn
            monomial_mask = tuple([int(var_name.Xn) for var_name in target_variables])
            support[monomial_mask] = support.get(monomial_mask, 0) + 1
        superpoly = []
        monomial_ntrails = dict()
        for monomial_mask in support.keys():            
            monomial = []
            for i in range(0, len(target_variables)):
                if monomial_mask[i] == 1:
                    monomial += [target_variables[i].varName]            
            monomial.sort()
            monomial = "*".join(monomial)
            if monomial == "":
                monomial = "1"
            monomial_ntrails[monomial] = support[monomial_mask]
            if support[monomial_mask] % 2 == 1:
                superpoly += [monomial]
        superpoly = " + ".join(superpoly)
        with open(self.result_file_name, "w") as outputfile:
            outputfile.write("Time used to solve the MILP problem: {:.02f}\n".format(elapsed_time))
            outputfile.write(f"Number of solutions: {number_of_solutions}\n") 
            for k, v in monomial_ntrails.items():
                outputfile.write("{:100}: {}\n".format(k, v))
            outputfile.write(superpoly)
        return superpoly

def parse_args():
    """
    parse input parameters
    """

    parser = ArgumentParser(description="This tool derives and solves the MILP "
                                        "model corresponding to integral analysis "
                                        "based on monomial prediction",
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-nr", "--nrounds", default=21, type=int, help="number of rounds\n")
    parser.add_argument("-tl", "--timelimit", default=5, type=int, help="set a time limit for the solver in seconds\n")
    parser.add_argument("-sp", "--superpoly", default=False, type=bool, help="extract the superpoly")
    parser.add_argument("-at", "--activeindices", nargs="+", default=[2], type=int, help="list of fixed input bits")
    parser.add_argument("-tr", "--targetbit", default=22, type=int, help="target output bit")
    return vars(parser.parse_args())

if __name__ == '__main__':
    locals().update(parse_args())
    warp = Warp(nrounds=nrounds)
    with open(warp.result_file_name, "w") as outputfile:
        outputfile.write(f"Results of applying monomial prediction method on {warp.nrounds} rounds of WARP\n")
    if superpoly:
        superpoly = warp.superpoly_recovery(active_indices=activeindices, target_output_bit=targetbit)
        print(superpoly)
    else:
        nb_dict = dict()
        for bit in range(0, 128):
            separator = "#"*100 + "\n"
            print(separator)
            print(f"Fixed index: {bit}")
            number_of_balanced_bits = warp.check_cube_method2(fixed_indices=[bit], time_limit=timelimit)
            nb_dict[bit] = number_of_balanced_bits
        print(nb_dict)
    
