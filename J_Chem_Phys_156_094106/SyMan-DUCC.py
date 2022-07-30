#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Notes: Any operators specific change to the code is marked with <CHANGE ME>
#        So, if you want to add define a new operator, you must check, not necessarily change,
#        all parts marked with <CHANGE ME>


# In[2]:


from sympy.combinatorics.partitions import Partition
from sympy.combinatorics.permutations import Permutation


# In[3]:


# The following declares the set of operators.
# Final form of the ordered operators is
# [Bra][Left][exp(-O)][observable][exp(O)][Right][Ket]

### Bra ###
Bra = ["E_V "]

### Ket ###
Ket = ["E^V "]

### Left ###
Left = [""]

### Right ###
Right = [""]

### Observable ###
Obs = ["F "]
# Obs = ["F ","V "]
# Obs = ["V "]

### Exponential ###
# Opp is the expression in the exponent written as a list.
# More specifically it is a list of 'expressions written as a list'
# Also, only the opperator for exp(O) is defined as exp(-O) will be taken care of internally
# Example 1: exp(T1-T1+) = ["T1 "," -T1+ "]
# Example 2: exp(T-T+), where T=T1+T2 = ["T1 ","-T1+ ","T2 ","-T2+ "] or ["T1 ","T2 ","-T1+ ","-T2+ "]

Opp = ["T1 ","-T1+ ","T2 ","-T2+ "]
exp_order=4 # Sets max number of operators per exponential, but also the commutator limit.
exp1 = [""]
exp2 = [""]

for i in range(0,exp_order):
    expanded_new=[""]
    for term in exp1:
        for o in Opp:
            # Flip the sign of operators since it is exp(-O) 
            if (o.count("-") == 1):
                o = o.replace('-','')
            else:
                o = "-"+o
            #
            if ((term+o).count("-") % 2) == 0:
                sign=""
            else:
                sign="-"
            expanded_new.append(sign+term.replace('-','')+o.replace('-',''))
    exp1=expanded_new[:]

for i in range(0,exp_order):
    expanded_new=[""]
    for term in exp2:
        for o in Opp:
            if ((term+o).count("-") % 2) == 0:
                sign=""
            else:
                sign="-"
            expanded_new.append(sign+term.replace('-','')+o.replace('-',''))
    exp2=expanded_new[:]


# In[4]:


# Form all combinations of 
operators_list=[]
for b in Bra:
    for l in Left:
        for e1 in exp1:
            for o in Obs:
                for e2 in exp2:
                    if((e1+e2).count("T")<=exp_order): #Limits the expansion
#                     if((e1+e2).count("T")==exp_order): #Limits the expansion
                        for r in Right:
                            for k in Ket:
                                if ((e1+e2).count("-") % 2) == 0:
                                    sign=""
                                else:
                                    sign="-"
                                operators_list.append(sign+" "+b+l+e1.replace('-','')+o+e2.replace('-','')+r+k)

# Uncomment for the full list of operator combinations before pruning
print(operators_list)


# In[5]:


from prettytable import PrettyTable
x = PrettyTable()
y = PrettyTable()
x.field_names = ["Expression", "# O+", "# V+", "# G+", "# O", "# V", "# G","Order"]
y.field_names = ["Expression", "# O+", "# V+", "# G+", "# O", "# V", "# G","Order","Error"]

evaluated_expressions=[]
for expression in operators_list:
    o_creation=0
    v_creation=0
    g_creation=0
    o_annihilation=0
    v_annihilation=0
    g_annihilation=0
    pt_order=0
    error=""

    expression_set= expression.split()

#   <CHANGE ME>
    for operator in expression_set:

        if(operator == "E_V"):
            v_annihilation += 1

        if(operator == "E^O"):
            o_creation += 1
            
        if(operator == "E^O_V"):
            o_creation += 1
            v_annihilation += 1
            
        if(operator == "E^OO_VV"):
            o_creation += 2
            v_annihilation += 2
            
        if(operator == "E^OO_V"):
            o_creation += 2
            v_annihilation += 1
            
        if(operator == "E^O_VV"):
            o_creation += 1
            v_annihilation += 2
            
        if(operator == "E^OO"):
            o_creation += 2
            
        if(operator == "E_VV"):
            v_annihilation += 2

        if(operator == "E^V"):
            v_creation += 1
            
        if(operator == "E_O"):
            o_annihilation += 1

        if(operator == "E_OO"):
            o_annihilation += 2
            
        if(operator == "E^VV"):
            v_creation += 2
            
        if(operator == "E^V_O"):
            v_creation += 1
            o_annihilation += 1
        
        if(operator == "F"):
            g_creation += 1
            g_annihilation += 1
            
        if(operator == "V"):
            g_creation += 2
            g_annihilation += 2
            pt_order += 1
            
        if(operator == "T1"):
            o_annihilation += 1
            v_creation += 1
            pt_order += 2

        if(operator == "T2"):
            o_annihilation += 2
            v_creation += 2
            pt_order += 1

        if(operator == "T1+"):
            v_annihilation += 1
            o_creation += 1
            pt_order += 2
            
        if(operator == "T2+"):
            v_annihilation += 2
            o_creation += 2
            pt_order += 1

    # Error Messages
    #   Check incoming/outgoing, left-extending/right-extending counts and remove those that
    #   cannot be fully contracted.
    if(o_creation > o_annihilation + g_annihilation):
        error+="*O+ > V + G*"
    if(v_creation > v_annihilation + g_annihilation):
        error+="*V+ > V + G*"
    if(o_annihilation > o_creation + g_creation):
        error+="*O > O+ + G+*"
    if(v_annihilation > v_creation + g_creation):
        error+="*V > V+ + G+*"

#   <CHANGE ME>
    # Remove combination if the first term kills the bra
#     if("-" in expression_set[0] and "E" in expression_set[1]):
#         if(expression_set[2] in ["T1","T2","T3","T4"]):
#             error+="*Second operator contains V+ & O*"
#     elif("-" in expression_set[0] and "E" not in expression_set[1]):
#         if(expression_set[1] in ["T1","T2","T3","T4"]):
#             error+="*First operator contains V+ & O*"
#     elif("-" not in expression_set[0] and "E" in expression_set[0]):
#         if(expression_set[1] in ["T1","T2","T3","T4"]):
#             error+="*Second operator contains V+ & O*"
#     elif("-" not in expression_set[0] and "E" not in expression_set[0]):
#         if(expression_set[0] in ["T1","T2","T3","T4"]):
#             error+="*First operator contains V+ & O*"

    # Remove combination if the last term kills the ket
    if("E" in expression_set[len(expression_set)-1]):
        if(expression_set[len(expression_set)-2] in ["T1+","T2+","T3+","T4+"]):
            error+="*2nd to Last operator contains V & O+*"
    elif(expression_set[len(expression_set)-1] in ["T1+","T2+","T3+","T4+"]):
        error+="*Last operator contains V & O+*"

    # Remove combination remove combination if 
#     if(expression.count("T")>exp_order):
#         error+="*Belongs to commutator > than exp_order*"  #outdated with added line in previous block
#
#   DONE FOR DUCC
#     if("T1+ T1 " in expression):
#         error+="*T1+ T1 in expression*"
#     if("T1+ T2 " in expression):
#         error+="*T1+ T2 in expression*"
#     if("T2+ T1 " in expression):
#         error+="*T2+ T1 in expression*"
#     if("T2+ T2 " in expression):
#         error+="*T2+ T2 in expression*"

    if(error == ""):        
        x.add_row([expression, o_creation, v_creation, g_creation, o_annihilation, v_annihilation, g_annihilation, pt_order])
        evaluated_expressions.append(expression) 
    else:
        y.add_row([expression, o_creation, v_creation, g_creation, o_annihilation, v_annihilation, g_annihilation, pt_order, error])

        
# print(evaluated_expressions)
print(x)
print(y)


# In[10]:


import math
import numpy as np
from ncon import ncon
from string import ascii_letters, digits
from fractions import Fraction
from IPython.core.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

def all_pairs(lst):
    if len(lst) < 2:
        yield []
        return
    if len(lst) % 2 == 1:
        # Handle odd length list
        for i in range(len(lst)):
            for result in all_pairs(lst[:i] + lst[i+1:]):
                yield result
    else:
        a = lst[0]
        for i in range(1,len(lst)):
            allow2 = True
            pair = (a,lst[i])

            if(identity_list[a]==identity_list[lst[i]]): allow2 = False
            if(("+" not in c_a_list[a]) and ("+" not in c_a_list[lst[i]])): allow2 = False
            if(("+" in c_a_list[a]) and ("+" in c_a_list[lst[i]])): allow2 = False
            if(identity_list[a] in is_excitation and identity_list[lst[i]] in is_excitation): allow2 = False
            if(identity_list[a] in is_deexcitation and identity_list[lst[i]] in is_deexcitation): allow2 = False
            if(("+" not in c_a_list[a]) and ("+" in c_a_list[lst[i]]) and ("o" in c_a_list[a])): allow2 = False
            if(("+" not in c_a_list[a]) and ("+" in c_a_list[lst[i]]) and ("o" in c_a_list[lst[i]])): allow2 = False
            if(("+" in c_a_list[a]) and ("+" not in c_a_list[lst[i]]) and ("v" in c_a_list[a])): allow2 = False
            if(("+" in c_a_list[a]) and ("+" not in c_a_list[lst[i]]) and ("v" in c_a_list[lst[i]])): allow2 = False

            if(allow2 == True):
                for rest in all_pairs(lst[1:i]+lst[i+1:]):
                    yield [pair] + rest


# Original
#         a = lst[0]
#         for i in range(1,len(lst)):
#             pair = (a,lst[i])
#             for rest in all_pairs(lst[1:i]+lst[i+1:]):
#                 yield [pair] + rest
            
                
def containsAny(str, set):
    for c in set:
        if c in str: return 1
    return 0

def compare_alphanumeric(first, second):
    for character in first:
        if character in ascii_letters + digits and character not in second:
            return False
    return True
                
############ Decalration section
occ_list=["i","j","k","l","m","n","o","p","q","r","M","N","O","P","Q","R"]
virt_list=["a","b","c","d","e","f","g","h","s","t","E","F","G","H","S","T"]
free_list=["e","f","g","h","m","n","o","p","q","r","M","N","O","P","Q","R","s","t","E","F","G","H","S","T"]
free_occ_list=["m","n","o","p","q","r","M","N","O","P","Q","R"]
free_virt_list=["e","f","g","h","s","t","E","F","G","H","S","T"]
fixed_occ_list=["i","j","k","l"]
fixed_virt_list=["a","b","c","d"]

nocc = 3
nvirt = 4
n = nocc + nvirt

F=np.random.rand(n,n)
V=np.random.rand(n,n,n,n)

EV=np.ones((nvirt))
EOdV=np.ones((nocc, nvirt))
EOdVOdV=np.random.rand(nocc, nvirt, nocc, nvirt)
EOdVOd=np.random.rand(nocc, nvirt, nocc)
EVOdV=np.random.rand(nvirt, nocc, nvirt)
EVV=np.random.rand(nvirt, nvirt)
EOd=np.ones((nocc))
EOdOd=np.random.rand(nocc, nocc)
EVd=np.ones((nvirt))
EO=np.ones((nocc))
EOO=np.random.rand(nocc, nocc)
EVdVd=np.random.rand(nvirt, nvirt)
EVdO=np.random.rand(nvirt, nocc)

T1=np.random.rand(nvirt, nocc)
T2=np.random.rand(nvirt, nocc, nvirt, nocc)

T1d=np.random.rand(nocc, nvirt)
T2d=np.random.rand(nocc, nvirt, nocc, nvirt)

for p in range(0,n):
    for q in range(0,n):
        for r in range (0, n):
            for s in range(0, n):
                if(p==r or q==s):
                    V[p,q,r,s] = 0.
                if(p<r and q<s):
                    element = V[p,q,r,s]
                    V[r,q,p,s] = -1.*element
                    V[p,s,r,q] = -1.*element
                    V[r,s,p,q] =  element

for i in range(0,nocc):
    for j in range(0,nocc):
        for a in range (0, nvirt):
            for b in range(0,nvirt):
                if(i==j or a==b):
                    T2[a,i,b,j] = 0.

                if(i<j and a<b):
                    element = T2[a,i,b,j]
                    T2[b,i,a,j] = -1.*element
                    T2[a,j,b,i] = -1.*element
                    T2[b,j,a,i] =  element      

for i in range(0,nocc):
    for j in range(0,nocc):
        for a in range (0, nvirt):
            for b in range(0,nvirt):
                if(i==j or a==b):
                    T2d[i,a,j,b] = 0.

                if(i<j and a<b):
                    element = T2d[i,a,j,b]
                    T2d[i,b,j,a] = -1.*element
                    T2d[j,a,i,b] = -1.*element
                    T2d[j,b,i,a] =  element       
                    
############
FullTable = PrettyTable()
FullTable.field_names = ["Expression", "Terms", "1/Weight", "Comm.", "PT Order", "HF or NHF", "Value"]
FullTable.align["Value"] = "r"
FullTable_Short=[]
for expression in evaluated_expressions:
    operator_list = expression.split()
    print(operator_list)
    
    c_a_string = ""
    has_bra = False
    has_ket = False
    is_excitation=[]
    is_deexcitation=[]
    # Identity is just a number assigned to an operator to make sure no contractions are within the operator
    if ("-" in operator_list[0]):
        identity = 1
    else:
        identity = 0
    identity_list = [] # Which operator a creation or annihilation operator belongs
    for operator in operator_list:
       
        if(operator == "E_V"):
            c_a_string += "v "
            identity_list.append(identity)
            identity += 1
            has_bra = True
            
        if(operator == "E^O_V"):
            c_a_string += "o+ v "
            identity_list.append(identity)
            identity_list.append(identity)
            identity += 1
            has_bra = True
            
        if(operator == "E^OO_VV"):
            c_a_string += "o+ v o+ v "
            identity_list.append(identity)
            identity_list.append(identity)
            identity_list.append(identity)
            identity_list.append(identity)
            identity += 1
            has_bra = True
            
        if(operator == "E^OO_V"):
            c_a_string += "o+ v o+ "
            identity_list.append(identity)
            identity_list.append(identity)
            identity_list.append(identity)
            identity += 1
            has_bra = True
            
        if(operator == "E^O_VV"):
            c_a_string += "v o+ v "
            identity_list.append(identity)
            identity_list.append(identity)
            identity_list.append(identity)
            identity += 1
            has_bra = True
            
        if(operator == "E^O"):
            c_a_string += "o+ "
            identity_list.append(identity)
            identity += 1
            has_bra = True
            
        if(operator == "E^OO"):
            c_a_string += "o+ o+ "
            identity_list.append(identity)
            identity_list.append(identity)
            identity += 1
            has_bra = True
            
        if(operator == "E_VV"):
            c_a_string += "v v "
            identity_list.append(identity)
            identity_list.append(identity)
            identity += 1
            has_bra = True

        if(operator == "E^V"):
            c_a_string += "v+ "
            identity_list.append(identity)
            identity += 1
            has_ket = True
            
        if(operator == "E_O"):
            c_a_string += "o "
            identity_list.append(identity)
            identity += 1
            has_ket = True
            
        if(operator == "E_OO"):
            c_a_string += "o o "
            identity_list.append(identity)
            identity_list.append(identity)
            identity += 1
            has_ket = True
            
        if(operator == "E^VV"):
            c_a_string += "v+ v+ "
            identity_list.append(identity)
            identity_list.append(identity)
            identity += 1
            has_ket = True
            
        if(operator == "E^V_O"):
            c_a_string += "v+ o "
            identity_list.append(identity)
            identity_list.append(identity)
            identity += 1
            has_ket = True
        
        if(operator == "F"):
            c_a_string += "g+ g "
            identity_list.append(identity)
            identity_list.append(identity)
            identity += 1
            
        if(operator == "V"):
            c_a_string += "g+ g g+ g "
            identity_list.append(identity)
            identity_list.append(identity)
            identity_list.append(identity)
            identity_list.append(identity)
            identity += 1
            
        if(operator == "T1"):
            c_a_string += "v+ o "
            is_excitation.append(identity)
            identity_list.append(identity)
            identity_list.append(identity)
            identity += 1

        if(operator == "T2"):
            c_a_string += "v+ o v+ o "
            is_excitation.append(identity)
            identity_list.append(identity)
            identity_list.append(identity)
            identity_list.append(identity)
            identity_list.append(identity)
            identity += 1

        if(operator == "T1+"):
            c_a_string += "o+ v "
            is_deexcitation.append(identity)
            identity_list.append(identity)
            identity_list.append(identity)
            identity += 1
            
        if(operator == "T2+"):
            c_a_string += "o+ v o+ v "
            is_deexcitation.append(identity)
            identity_list.append(identity)
            identity_list.append(identity)
            identity_list.append(identity)
            identity_list.append(identity)
            identity += 1

    c_a_list=c_a_string.split()
    print(c_a_string)
#     print(identity_list)
        
    length_list = range(0,len(c_a_list))
    allowed_pairs_list = []
    for x in all_pairs(list(length_list)):
        allow = True
#         print(x)
        # Check all contractions
        for contracted_pair in x:
#             # Remove contraction if between two excitation operators
#             if(identity_list[contracted_pair[0]] in is_excitation and identity_list[contracted_pair[1]] in is_excitation):
#                 allow= False
#             if not allow: break

#             # Remove contraction if between two excitation operators
#             if(identity_list[contracted_pair[0]] in is_deexcitation and identity_list[contracted_pair[1]] in is_deexcitation):
#                 allow= False
#             if not allow: break
                        
#             # Remove combination if there is a contraction within an operator
#             if(identity_list[contracted_pair[0]] == identity_list[contracted_pair[1]]):
# #                 print("contraction withing same operator")
#                 allow= False
#             if not allow: break
#
#             # Remove combination if contraction is between two annihilation
#             if(("+" not in c_a_list[contracted_pair[0]]) and ("+" not in c_a_list[contracted_pair[1]])):
# #                 print("contraction between two annihilation")
#                 allow= False
#             if not allow: break
# #
#             # Remove combination if contraction is between two creation
#             if(("+" in c_a_list[contracted_pair[0]]) and ("+" in c_a_list[contracted_pair[1]])):
# #                 print("contraction between two creation")
#                 allow= False
#             if not allow: break
#
#             # x-x+ contraction con only be between virtual orbitals
#             if(("+" not in c_a_list[contracted_pair[0]]) and ("+" in c_a_list[contracted_pair[1]]) and ("o" in c_a_list[contracted_pair[0]])):
# #                 print("x-x+ where x is occ.")
#                 allow= False
#             if not allow: break
#             if(("+" not in c_a_list[contracted_pair[0]]) and ("+" in c_a_list[contracted_pair[1]]) and ("o" in c_a_list[contracted_pair[1]])):
# #                 print("x-x+ where x+ is occ.")
#                 allow= False
#             if not allow: break
# #
#             # x-x+ contraction con only be between occupied orbitals
#             if(("+" in c_a_list[contracted_pair[0]]) and ("+" not in c_a_list[contracted_pair[1]]) and ("v" in c_a_list[contracted_pair[0]])):
# #                 print("x+-x where x+ is virt.")
#                 allow= False
#             if not allow: break
#             if(("+" in c_a_list[contracted_pair[0]]) and ("+" not in c_a_list[contracted_pair[1]]) and ("v" in c_a_list[contracted_pair[1]])):
# #                 print("x+-x where x is virt.")
#                 allow= False
#             if not allow: break
#            
            if(has_bra and has_ket):
                if("-" in operator_list[0]):
                    if((identity_list[contracted_pair[0]] == 1) and (identity_list[contracted_pair[1]] == identity-1)):
#                         print("Contraction between bra and ket")
                        allow= False
                        if not allow: break
                elif((identity_list[contracted_pair[0]] == 0) and (identity_list[contracted_pair[1]] == identity-1)):
#                         print("Contraction between bra and ket")
                        allow= False
                        if not allow: break
#
        if allow:
            allowed_pairs_list.append(x)
            
    print("******************")
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
    ExpressionTable = PrettyTable()
    ExpressionTable.field_names = ["Terms", "1/Weight", "PT Order", "HF or NHF", "Value"]
    ExpressionTable.align["Value"] = "r"
    
    for allowed_contraction_set in allowed_pairs_list:
        
#         print(allowed_contraction_set)
        free_occ_index= 0
        free_virt_index= 0
        fixed_occ_index= 0
        fixed_virt_index= 0
        
        # Compute the Parity of the permutation
        permutation_list=[]
        for contracted_pair in allowed_contraction_set:
            permutation_list.append(contracted_pair[0])
            permutation_list.append(contracted_pair[1])
#         print(permutation_list)
        p = Permutation(permutation_list)
        
        # Final sign depnds on the sign of the expression and 
        if((p.parity() == 1)  and ("-" in operator_list[0])):
            sign = 1
        if((p.parity() == 1)  and ("-" not in operator_list[0])):
            sign = -1
        if((p.parity() == 0)  and ("-" in operator_list[0])):
            sign = -1
        if((p.parity() == 0)  and ("-" not in operator_list[0])):
            sign = 1
#         print("Sign = ", sign)
        
        
        bra_index = -1
        ket_index = -1
        if(has_bra):
            if ("-" in operator_list[0]):
                bra_index = 1
            else:
                bra_index = 0
        if(has_ket):
                ket_index = len(operator_list)-1
        
        labeled_string = "X" * len(c_a_list)
        for contracted_pair in allowed_contraction_set:
            if((identity_list[contracted_pair[0]] == bra_index) and ("v" in c_a_list[contracted_pair[0]])):
                labeled_string = labeled_string[:contracted_pair[0]] + fixed_virt_list[fixed_virt_index] + labeled_string[contracted_pair[0]+1:]
                labeled_string = labeled_string[:contracted_pair[1]] + fixed_virt_list[fixed_virt_index] + labeled_string[contracted_pair[1]+1:]
                fixed_virt_index += 1
            elif((identity_list[contracted_pair[0]] == bra_index) and ("o" in c_a_list[contracted_pair[0]])):
                labeled_string = labeled_string[:contracted_pair[0]] + fixed_occ_list[fixed_occ_index] + labeled_string[contracted_pair[0]+1:]
                labeled_string = labeled_string[:contracted_pair[1]] + fixed_occ_list[fixed_occ_index] + labeled_string[contracted_pair[1]+1:]
                fixed_occ_index += 1
            elif((identity_list[contracted_pair[1]] == ket_index) and ("v" in c_a_list[contracted_pair[1]])):
                labeled_string = labeled_string[:contracted_pair[0]] + fixed_virt_list[fixed_virt_index] + labeled_string[contracted_pair[0]+1:]
                labeled_string = labeled_string[:contracted_pair[1]] + fixed_virt_list[fixed_virt_index] + labeled_string[contracted_pair[1]+1:]
                fixed_virt_index += 1
            elif((identity_list[contracted_pair[1]] == ket_index) and ("o" in c_a_list[contracted_pair[1]])):
                labeled_string = labeled_string[:contracted_pair[0]] + fixed_occ_list[fixed_occ_index] + labeled_string[contracted_pair[0]+1:]
                labeled_string = labeled_string[:contracted_pair[1]] + fixed_occ_list[fixed_occ_index] + labeled_string[contracted_pair[1]+1:]
                fixed_occ_index += 1
            elif(("v" in c_a_list[contracted_pair[0]]) or ("v" in c_a_list[contracted_pair[1]])):
                labeled_string = labeled_string[:contracted_pair[0]] + free_virt_list[free_virt_index] + labeled_string[contracted_pair[0]+1:]
                labeled_string = labeled_string[:contracted_pair[1]] + free_virt_list[free_virt_index] + labeled_string[contracted_pair[1]+1:]
                free_virt_index += 1
            elif(("o" in c_a_list[contracted_pair[0]]) or ("o" in c_a_list[contracted_pair[1]])):
                labeled_string = labeled_string[:contracted_pair[0]] + free_occ_list[free_occ_index] + labeled_string[contracted_pair[0]+1:]
                labeled_string = labeled_string[:contracted_pair[1]] + free_occ_list[free_occ_index] + labeled_string[contracted_pair[1]+1:]
                free_occ_index += 1

        # put brackets around labels for given operators
        labeled_string = "{" + labeled_string
        position = 2
        for identity in range(0,len(identity_list)-1):
            if(identity_list[identity] != identity_list[identity+1]):
                labeled_string = labeled_string[:position] + "}  {" + labeled_string[position:]
                position += 4
            position +=1
                
        labeled_string = " " + labeled_string + "}"
        
        start = 0
        count = 0
        if("-" in operator_list[0]):
            start = 1
            
        for operator in operator_list[start:]:
            labeled_string = labeled_string[:labeled_string.index(" {")] + operator + "{" + labeled_string[labeled_string.index(" {")+2:]
            
#         print(labeled_string)

###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
        
        split_labeled_string = labeled_string.split()
        TensorArray = []
        IndexArray = []

# *** Depending on the method, this may need to be commented out
# SKIP EVALUATION IF THERE IS SOME DISCONNECTED PART
        # Skip evaluation if F or V is disconnected.
        skip_eval = False
        if("T" in labeled_string):
            for term in split_labeled_string:
                if("F{" in term and containsAny(term,free_list)==False):
#                     print("skipping evaluation as F is disconnected")
                    skip_eval = True
                elif("V{" in term and "E" not in term and containsAny(term,free_list)==False):
#                     print("skipping evaluation as V is disconnected")
                    skip_eval = True

#### This part can be made more general.
#### With up to 3 commutators, disconnected terms involving two or three operators
        #Look for disconnected term involving two elements
        if(skip_eval == False):
            for term_a in split_labeled_string:
                str1 = "y"
                if("E" not in term_a):
                    str1 = term_a[term_a.index("{")+1:term_a.index("}")]
                for term_b in split_labeled_string[split_labeled_string.index(term_a)+1:]:
                    str2 = "z"
                    if("E" not in term_b):
                        str2 = term_b[term_b.index("{")+1:term_b.index("}")]
                    if(compare_alphanumeric(str1, str2) == True and compare_alphanumeric(str2, str1) == True):
#                         print(labeled_string," has a disconnected term")
                        skip_eval = True
        # Disconnected term involving 3 operators is the same as 1 operator having only fixed indices. With the exception of having only F or V
        if((skip_eval == False) and ("T" in labeled_string)):
            for term_a in split_labeled_string:
                has_free = False
                if("E" not in term_a):
                    for index in free_list:
                        if(index in term_a[term_a.index("{")+1:term_a.index("}")]):
                            has_free = True
                    if(has_free == False):
                        skip_eval = True
#                         print("Skipping because of ",term_a)
    
    
        if(skip_eval == False):
            weight = 1
            pt_order=0
            NHF_Term = ""
            for term in split_labeled_string:
                TermLabelArray = []
                
                if("E_V{" in term):
                    TensorArray.append(EV)
                if("E^O_V{" in term):
                    TensorArray.append(EOdV)
                if("E^OO_VV{" in term):
                    TensorArray.append(EOdVOdV)
                if("E^OO_V{" in term):
                    TensorArray.append(EOdVOd)
                if("E^O_VV{" in term):
                    TensorArray.append(EVOdV)
                if("E^O{" in term):
                    TensorArray.append(EOd)
                if("E^OO{" in term):
                    TensorArray.append(EOdOd)
                if("E_VV{" in term):
                    TensorArray.append(EVV)
                if("E^V{" in term):
                    TensorArray.append(EVd)
                if("E_O{" in term):
                    TensorArray.append(EO)
                if("E_OO{" in term):
                    TensorArray.append(EOO)
                if("E^VV{" in term):
                    TensorArray.append(EVdVd)
                if("E^V_O{" in term):
                    TensorArray.append(EVdO)
                if("F{" in term):
                    if(term[2] in occ_list and term[3] in occ_list):
                        TensorArray.append(F[0:nocc,0:nocc])
                    elif(term[2] in virt_list and term[3] in occ_list):
                        TensorArray.append(F[nocc:n,0:nocc])
                    elif(term[2] in occ_list and term[3] in virt_list):
                        TensorArray.append(F[0:nocc,nocc:n])
                    elif(term[2] in virt_list and term[3] in virt_list):
                        TensorArray.append(F[nocc:n,nocc:n])
                    if(term[term.index("{")+1] in occ_list and term[term.index("{")+2] in virt_list):
                        NHF_Term = "Non-HF"
                    if(term[term.index("{")+1] in virt_list and term[term.index("{")+2] in occ_list):
                        NHF_Term = "Non-HF"
                if("V{" in term and "E" not in term):
                    weight *= 4
                    pt_order += 1
                    if(term[2] in occ_list and term[3] in occ_list and term[4] in occ_list and term[5] in occ_list):
                        TensorArray.append(V[0:nocc,0:nocc,0:nocc,0:nocc])
                    elif(term[2] in occ_list and term[3] in occ_list and term[4] in occ_list and term[5] in virt_list):
                        TensorArray.append(V[0:nocc,0:nocc,0:nocc,nocc:n])
                    elif(term[2] in occ_list and term[3] in occ_list and term[4] in virt_list and term[5] in occ_list):
                        TensorArray.append(V[0:nocc,0:nocc,nocc:n,0:nocc])
                    elif(term[2] in occ_list and term[3] in occ_list and term[4] in virt_list and term[5] in virt_list):
                        TensorArray.append(V[0:nocc,0:nocc,nocc:n,nocc:n])
                    elif(term[2] in occ_list and term[3] in virt_list and term[4] in occ_list and term[5] in occ_list):
                        TensorArray.append(V[0:nocc,nocc:n,0:nocc,0:nocc])
                    elif(term[2] in occ_list and term[3] in virt_list and term[4] in occ_list and term[5] in virt_list):
                        TensorArray.append(V[0:nocc,nocc:n,0:nocc,nocc:n])
                    elif(term[2] in occ_list and term[3] in virt_list and term[4] in virt_list and term[5] in occ_list):
                        TensorArray.append(V[0:nocc,nocc:n,nocc:n,0:nocc])
                    elif(term[2] in occ_list and term[3] in virt_list and term[4] in virt_list and term[5] in virt_list):
                        TensorArray.append(V[0:nocc,nocc:n,nocc:n,nocc:n])
                    elif(term[2] in virt_list and term[3] in occ_list and term[4] in occ_list and term[5] in occ_list):
                        TensorArray.append(V[nocc:n,0:nocc,0:nocc,0:nocc])
                    elif(term[2] in virt_list and term[3] in occ_list and term[4] in occ_list and term[5] in virt_list):
                        TensorArray.append(V[nocc:n,0:nocc,0:nocc,nocc:n])
                    elif(term[2] in virt_list and term[3] in occ_list and term[4] in virt_list and term[5] in occ_list):
                        TensorArray.append(V[nocc:n,0:nocc,nocc:n,0:nocc])
                    elif(term[2] in virt_list and term[3] in occ_list and term[4] in virt_list and term[5] in virt_list):
                        TensorArray.append(V[nocc:n,0:nocc,nocc:n,nocc:n])
                    elif(term[2] in virt_list and term[3] in virt_list and term[4] in occ_list and term[5] in occ_list):
                        TensorArray.append(V[nocc:n,nocc:n,0:nocc,0:nocc])
                    elif(term[2] in virt_list and term[3] in virt_list and term[4] in occ_list and term[5] in virt_list):
                        TensorArray.append(V[nocc:n,nocc:n,0:nocc,nocc:n])
                    elif(term[2] in virt_list and term[3] in virt_list and term[4] in virt_list and term[5] in occ_list):
                        TensorArray.append(V[nocc:n,nocc:n,nocc:n,0:nocc])
                    elif(term[2] in virt_list and term[3] in virt_list and term[4] in virt_list and term[5] in virt_list):
                        TensorArray.append(V[nocc:n,nocc:n,nocc:n,nocc:n])
                if("T1{" in term):
                    pt_order += 2
                    TensorArray.append(T1)
                if("T2{" in term):
                    weight *= 4
                    pt_order += 1
                    TensorArray.append(T2)
                if("T1+{" in term):
                    pt_order += 2
                    TensorArray.append(T1d)
                if("T2+{" in term):
                    weight *= 4
                    pt_order += 1
                    TensorArray.append(T2d)
            
                for label in range(term.index("{")+1,term.index("}")):
                    TermLabelArray.append(ord(term[label])-96)
                
                IndexArray.append(TermLabelArray)
            
            ### Note about effective Hamiltonians: When doing wick's theorem on string of operators for an effective Hamiltonian,
            ### You have to remember that the remaining fixed lines will be reordered such that Y+ operators are before Y. Without 
            ### coding this up, you may have to change the sign with the following line:
#             sign *= -1

            ### In some cases, the the indicies for the bra or ket need to be swapped in the whole string (if you did only the bra
            ### or ket there would be a sign change). 
#             if("E^O_VV{aib}" in labeled_string):
#                 labeled_string = labeled_string.replace("a", "z")
#                 labeled_string = labeled_string.replace("b", "a")
#                 labeled_string = labeled_string.replace("z", "b")
                

            E = sign*ncon(TensorArray,IndexArray)

            t_count = 0
            element_count = 0
            for term in split_labeled_string:
                if("T" in term):
                    t_count += 1
                    element_count += 1
                if("F" in term or ("V" in term and "E" not in term)):
                    weight *= math.factorial(t_count)
                    t_count = 0
                    element_count += 1
                if(term == split_labeled_string[len(split_labeled_string)-1]):
                    weight *= math.factorial(t_count)
                            
            if(sign==-1):
                ExpressionTable.add_row(["- "+labeled_string,weight,pt_order,NHF_Term,"{:.12f}".format(E)])
                FullTable.add_row([expression,"- "+labeled_string,weight,element_count-1,pt_order,NHF_Term,"{:.12f}".format(E)])
                FullTable_Short.append([expression,"- "+labeled_string,weight,element_count-1,pt_order,NHF_Term,"{:.10f}".format(E)])
            else:
                ExpressionTable.add_row(["  "+labeled_string,weight,pt_order,NHF_Term,"{:.12f}".format(E)])
                FullTable.add_row([expression,"  "+labeled_string,weight,element_count-1,pt_order,NHF_Term,"{:.12f}".format(E)])
                FullTable_Short.append([expression,"  "+labeled_string,weight,element_count-1,pt_order,NHF_Term,"{:.10f}".format(E)])

    print(ExpressionTable)
    

    print()
    print("-----------------------------------------------------------------------------------")
print(FullTable)


ShortFullTable = PrettyTable()
ShortFullTable.field_names = ["Expression", "Terms", "Weight", "Comm.", "PT Order", "HF or NHF", "Value"]
ShortFullTable.align["Value"] = "r"
ShortTableList=[]
for tabletuple1 in FullTable_Short:
    count = 1
    additionalterm_list = []
    for tabletuple2 in FullTable_Short[FullTable_Short.index(tabletuple1)+1:]:
        if(tabletuple1[6]==tabletuple2[6] and tabletuple1[0]==tabletuple2[0]):
            count += 1
            additionalterm_list.append(tabletuple2)
    for term in additionalterm_list:
        FullTable_Short.remove(term)
    ShortFullTable.add_row([tabletuple1[0],tabletuple1[1],1/(tabletuple1[2]/count),tabletuple1[3],tabletuple1[4],tabletuple1[5],tabletuple1[6]])
    ShortTableList.append([tabletuple1[0],tabletuple1[1],1/(tabletuple1[2]/count),tabletuple1[3],tabletuple1[4],tabletuple1[5],tabletuple1[6],abs(float(tabletuple1[6]))])
print(ShortFullTable)
ShortTableList.sort(key=lambda a:a[7])
ShortTableList.sort(key=lambda a:a[4])
ShortTableList.sort(key=lambda a:a[3])

ShortFullTable_sorted = PrettyTable()
ShortFullTable_sorted.field_names = ["Expression", "Terms", "Weight", "Comm.", "PT Order", "HF or NHF", "Value"]
ShortFullTable_sorted.align["Value"] = "r"
for tabletuple in ShortTableList:
    ShortFullTable_sorted.add_row([tabletuple[0],tabletuple[1],tabletuple[2],tabletuple[3],tabletuple[4],tabletuple[5],tabletuple[6]])
print(ShortFullTable_sorted)

FinalTable2 = PrettyTable()
FinalTable2.field_names = ["Expression", "Terms", "Weight", "Comm.", "PT Order", "HF or NHF", "Value"]
FinalTable2.align["Value"] = "r"
for i in range (0,exp_order+1):
    FinalTable = PrettyTable()
    FinalTable.field_names = ["Expression", "Terms", "Weight", "Comm.", "PT Order", "HF or NHF", "Value"]
    FinalTable.align["Value"] = "r"
    current_value = 0 
    for tabletuple1 in ShortTableList:
        if(int(tabletuple1[3])==i):
            if(float(tabletuple1[7])==current_value):
                continue
            count = 0
            current_value = float(tabletuple1[7])
            if(float(tabletuple1[6])<0):
                weight_sum = -1.0*float(tabletuple1[2])
            else:
                weight_sum = float(tabletuple1[2])
            for tabletuple2 in ShortTableList[ShortTableList.index(tabletuple1)+1:]:
                if(tabletuple1[7] == tabletuple2[7]):
                    count += 1
            for tabletuple2 in ShortTableList[ShortTableList.index(tabletuple1)+1:ShortTableList.index(tabletuple1)+1+count]:
                if(float(tabletuple2[6])<0):
                    weight_sum += -1.0*float(tabletuple2[2])
                else:
                    weight_sum += float(tabletuple2[2])
            isdone = False
            for tabletuple in ShortTableList[ShortTableList.index(tabletuple1):ShortTableList.index(tabletuple1)+1+count]:
                elements = tabletuple[0].split()
                if("F" in elements):
                    if("T1+" in elements[elements.index("F"):] or "T2+" in elements[elements.index("F"):]):
                        continue
                if("V" in elements):
                    if("T1+" in elements[elements.index("V"):] or "T2+" in elements[elements.index("V"):]):
                        continue
                if("F" in elements):
                    if("T1" in elements[:elements.index("F")] or "T2" in elements[:elements.index("F")]):
                        continue
                if("V" in elements):
                    if("T1" in elements[:elements.index("V")] or "T2" in elements[:elements.index("V")]):
                        continue
                if(isdone == False):
                    FinalTable.add_row([tabletuple[0],tabletuple[1],abs(weight_sum),tabletuple[3],tabletuple[4],tabletuple[5],tabletuple[6]])
                    FinalTable2.add_row([tabletuple[0],tabletuple[1],abs(weight_sum),tabletuple[3],tabletuple[4],tabletuple[5],tabletuple[6]])
                    isdone == True
    print(FinalTable)   
print(FinalTable2)


# In[ ]:





# In[ ]:





# In[ ]:




