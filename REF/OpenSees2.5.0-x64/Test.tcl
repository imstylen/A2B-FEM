wipe
# Create ModelBuilder (with two-dimensions and 3 DOF/node)
model basic -ndm 2 -ndf 3


# Create nodes
#    tag        X       Y 
node  1       0.0       0.0 
node  2       0.0       120.0 
node  3       120.0     120.0
node  4       120.0     0.0 
node  5       -120.0     0.0 

# Fix supports at base of columns
#    tag   DX   DY   RZ
fix   1     1    1    1
fix   4     1    1    1


#                tag 
geomTransf Linear 1  

# Create the beam elements
#                          tag ndI ndJ     A       E    Iz   transfTag
element elasticBeamColumn   1   1   2    3.54    29000  53.8    1
element elasticBeamColumn   2   2   3    3.54    29000  53.8    1
element elasticBeamColumn   3   3   4    3.54    29000  53.8    1
element elasticBeamColumn   4   1   5    3.54    29000  53.8    1


# Define gravity loads
# --------------------

# Create a Plain load pattern with a Linear TimeSeries
pattern Plain 1 "Linear" {

        # Create nodal loads at nodes 3 & 4
	#    nd    FX   FY      MZ 
	load  2   10    0.0     0.0
	load  5   0.0   -1.0     0
    
}

# ------------------------------
# End of model generation
# ------------------------------



# ------------------------------
# Start of analysis generation
# ------------------------------

# Create the system of equation, a sparse solver with partial pivoting
system FullGeneral

# Create the constraint handler, the transformation method
constraints Transformation

# Create the DOF numberer, the reverse Cuthill-McKee algorithm
numberer RCM

# Create the convergence test, the norm of the residual with a tolerance of 
# 1e-12 and a max number of iterations of 10
test NormDispIncr 1.0e-12  10 3

# Create the solution algorithm, a Newton-Raphson algorithm
algorithm Newton

# Create the integration scheme, the LoadControl scheme using steps of 0.1 
integrator LoadControl 1

# Create the analysis object
analysis Static

# ------------------------------
# End of analysis generation
# ------------------------------


# ------------------------------
# Finally perform the analysis
# ------------------------------

recorder Node -file nodesRxn.out -node 1 4 -dof 1 2 3 reaction;
recorder Element -file ElementGlobalForce.out -ele 1 2 3 4 globalForce;
recorder Element -file ElementLocalForce.out -ele 1 2 3 4 localForce;

# perform the gravity load analysis, requires 10 steps to reach the load level
analyze 1
printA
# printA -file "stiffness.out"

# Print out the state of nodes 3 and 4

# Print out the state of element 1
#print ele 1 
