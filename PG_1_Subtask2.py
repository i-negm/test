"""
This code is a modification of PG_1_Subtask1 intended to include a changing space charge
at the surface of the conductor as well as the Corona discharge voltage
"""
from xml.etree.ElementTree import QName
from dolfin import *
from fenics import *
from mshr import *
import scipy.constants as c
import numpy as np
import matplotlib.pyplot as plt



#
# Functions
#

def calculate_phi(MPhi, rhoS):
    epsi = Constant(c.epsilon_0)        
    checkIter=False
    n = 0
    while not checkIter :
        n = n+1
        print(n, "iteration of steps 2, 3, 4")
        #changed in 0.1.3.C
        if n > LoopBreakPoint :
            print('more then ', LoopBreakPoint, ' iterations')
            break
            # stops the loop if we run it to often

        #changed in 0.1.3.C
        phi_old = MPhi.compute_vertex_values()

        if True:
            phi_current = MPhi.compute_vertex_values()
            minVal = np.amin(phi_current)
            print("minimal Value = ",minVal)
            #input("Press Enter to continue...")

        # step 2 Generate gradient of Phi, should be a vector 
        E = Function(V)  
        E.assign(project(-grad(MPhi), V))
            
        print('step 2 finished')
        if Test:
            IterationE=plot(E, title='show E')
            plt.colorbar(IterationE)
            plt.show()

        #################################################################
        #calculating rho0
        marker = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)
        for k in facets(mesh):
            marker[k] = 0.5 - DOLFIN_EPS < k.midpoint().x() < 0.5 + DOLFIN_EPS

        ## Create submesh ##
        submesh = MeshView.create(marker, 1)        

        if True:


            LM = FunctionSpace(submesh, "Real", 0)
            W = MixedFunctionSpace(Q,LM)
        else:
            P = FiniteElement("Lagrange", mesh.ufl_cell(), 2)
            LM = FiniteElement("Real", mesh.ufl_cell(), 0)
            W = FunctionSpace(mesh, MixedElement(P, LM))

        Trial1=Function(W)
        rho0,_=Trial1.split()
        (s,s_lm) = TestFunctions(W)

        ##########################################################################
        # different ways to generate an Intitial Guess
        if n==1 and False:
            # initial guess through solving of a simplified PDE and DirichletBC
            #rho0 = Constant(10e-6)
            #needs to be formula since differentiation = 0 
            guess0 = TestFunction(Q)
            h0= TrialFunction(Q)

            bla = Constant(1)
            origin = 'near(x[0], 0.49)'
            oriFunc= Expression('50*x[1]*x[1]', degree=1)
            #bo = DirichletBC(Q, Constant(5),origin )
            bo = DirichletBC(Q, oriFunc,origin )


            guess0 = Function(Q)

            test0= -1*(guess0**2/epsi)*h0*ds(1)

            test0J=derivative(test0,guess0)
            print(test0J)

            test0problem = NonlinearVariationalProblem(test0, guess0, bo , test0J)
            test0solver = NonlinearVariationalSolver(test0problem)

            test0prm = test0solver.parameters

            test0solver.solve()

            print(guess0)

            
        elif n==1 and False:
            rho0=np.random.rand(1) #random intitial guess
            print('nothing')


        ################################################################################
        # Different ways to Calculate rho0
        if False:
            class Conductor(SubDomain):
                def inside(self, x, on_boundary):
                    #return on_boundary and x[0]>0.48 and x[0]<0.51 and x[1]>0.54 and x[1]<0.56
                    return (between(x[0], (electrodePosX-CondSize, electrodePosX+CondSize)) and between(x[1], (electrodePosY-CondSize, electrodePosY+CondSize)))


            conductor = Conductor()

            boundarieCon = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
            boundarieCon.set_all(0)
            conductor.mark(boundarieCon, 1)

            dBound = ds(subdomain_data=boundarieCon) 

            F0=(rho0**2/epsi)*s*dBound(1) + (voltage-MPhi)*s_lm*dBound(1)

            rhoJ = derivative(F0, rho0, s)

            rho0=Function(Q)

            rho0problem = NonlinearVariationalProblem(F0, rho0, J=rho0J )

            rho0solver = NonlinearVariationalSolver(rho0problem)

            rho0prm =rho0solver.parameters
            rho0prm['newton_solver']['absolute_tolerance'] = 1E-5
            rho0prm['newton_solver']['relative_tolerance'] = 1E-5
            rho0prm['newton_solver']['maximum_iterations'] = 5 
            rho0prm['newton_solver']['relaxation_parameter'] = 1.0
            rho0solver.solve()

            rho0_values=rho0.compute_vertex_values()
            rho0=rho0_values[minDisIndex]
            
        elif False:
            dSurf0 = Measure("ds", domain=W.sub_space(0).mesh())
            dSurf1 = Measure("ds", domain=W.sub_space(1).mesh())

            F0=(rho0**2/epsi)*s*dSurf0 + (voltage-MPhi)*s_lm*dSurf1

            rhoJ = derivative(F0, rho0, s)

            rho0problem = NonlinearVariationalProblem(F0, Trial1, None, rho0J )

            rho0J=derivative(F0,rho0,s)
            
            
            rho0=Function(Q)

            rho0problem = NonlinearVariationalProblem(F0, rho0, J=rho0J )

            rho0solver = NonlinearVariationalSolver(rho0problem)

            rho0prm =rho0solver.parameters
            rho0prm['newton_solver']['absolute_tolerance'] = 1E-5
            rho0prm['newton_solver']['relative_tolerance'] = 1E-5
            rho0prm['newton_solver']['maximum_iterations'] = 5 
            rho0prm['newton_solver']['relaxation_parameter'] = 1.0
            rho0solver.solve()

            rho0_values=rho0.compute_vertex_values()
            rho0=rho0_values[minDisIndex]

        elif False:

            dField0 = Measure("dx", domain=W.sub_space(0).mesh())
            dField1 = Measure("dx", domain=W.sub_space(1).mesh())

            F0=(rho0**2/epsi)*s*dField0 + (voltage-MPhi)*s_lm*dField1

            rhoJ = derivative(F0, rho0, s)

            rho0problem = NonlinearVariationalProblem(F0, Trial1, None, rho0J )

            rho0J=derivative(F0,rho0,s)
            
            
            rho0=Function(Q)

            rho0problem = NonlinearVariationalProblem(F0, rho0, J=rho0J )

            rho0solver = NonlinearVariationalSolver(rho0problem)

            rho0prm =rho0solver.parameters
            rho0prm['newton_solver']['absolute_tolerance'] = 1E-5
            rho0prm['newton_solver']['relative_tolerance'] = 1E-5
            rho0prm['newton_solver']['maximum_iterations'] = 5 
            rho0prm['newton_solver']['relaxation_parameter'] = 1.0
            rho0solver.solve()

            rho0_values=rho0.compute_vertex_values()
            rho0=rho0_values[minDisIndex]

        else:
            # Backup: set rho0 as constant
            rho0=Constant(1*10E-6)
            

    ###################################################################################

        if False:
            # other solvers tried

            F0= (rho0**2/epsi)*s*ds(1) + (voltage-MPhi)*s_lm*dx #+ Constant(0)*rho0*s*dx

            if False :
                params = {'nonlinear_solver': 'snes',
                'snes_solver':
                    {
                        'linear_solver'           : 'mumps',
                        'absolute_tolerance'      : 1e-10,
                        'relative_tolerance'      : 1e-10,
                        'maximum_iterations'      : 20,
                    }
                }
                rho0solver = NonlinearVariationalSolver(rho0problem)

                rho0solver.parameters.update(params)
            elif False:
                rho0solver = NonlinearVariationalSolver(rho0problem)

                rho0prm = rho0solver.parameters
                rho0prm['newton_solver']['absolute_tolerance'] = 1E-10
                rho0prm['newton_solver']['relative_tolerance'] = 1E-10
                rho0prm['newton_solver']['maximum_iterations'] = 5
                rho0prm['newton_solver']['relaxation_parameter'] = 1.0
        

            #rho0solver.solve()
            print("rho0 =",rho0)


        print ('rho0 calculated')
        #################################################################
        # calculating rho_Space
        if (rhoS == None):
            rhoS=TrialFunction(Q)
            rhoS=Function(Q)
        u=TestFunction(Q)


        h = Constant(0)
        i = Constant(0)

        if n==1 and False :
            initiaGuessFile = HDF5File(MPI.comm_world,"static.h5","r")
            initiaGuessFile.read(rhoS,"rho-static")
            initiaGuessFile.close()
            
            print('initial guess loaded')
        
        #F = (dot((rhoS**2) / epsi, u) - dot(dot(E, grad(rhoS)), u))*dx
        FS= dot(((rho0+rhoS)**2)/epsi - dot(E,grad(rhoS)),u)*dx
        LC = h*u*dx #+ i*u*ds(1)


        rhoJ=derivative(FS,rhoS)
        if Test:
            print('derivate calculated, starting solver')


        rhoSproblem = NonlinearVariationalProblem(FS, rhoS, bcd ,rhoJ)

        rhoSsolver = NonlinearVariationalSolver(rhoSproblem)

        rhoSprm = rhoSsolver.parameters

        rhoSsolver.solve()

        print('calculated rho')
        if Test:
            IterationRhoS=plot(rhoS, title='show rhoS')
            plt.colorbar(IterationRhoS)
            plt.show()


        ################################################################################
        # Calculating the Poisson Equation

        #step 4
        # Equation 9
        MPhi = TrialFunction(Q)
        w = TestFunction(Q)

        aP = inner(grad(MPhi), grad(w))*dx 
        LP = -((rho0+rhoS)**2/epsi)*w*dx - Ec*w*ds(1)

        MPhi = Function(Q)

        solve(aP == LP, MPhi,bcp_ground) 
        print('calculated MPhi')


        # Plot solution
        import matplotlib.pyplot as plt
        if Test:
            print('Result MPhi in itteration ', n)

            IterationPhi=plot(MPhi, title='Result MPhi in itteration')
            plt.colorbar(IterationPhi)
            plt.show()
            

        # generating new checkIter
        phi_current = MPhi.compute_vertex_values()
        checkIter = (abs(phi_current - phi_old)< tol ).all()
        # should compare the values of the old and new MPhi at the vertices 
        print('finished itteration ', n)


        if n < minIter and checkIter == True:
            checkIter = False
            print('loop has ended before minimum ammount of iterations needed')
        if n >= maxIter and checkIter == False:
            checkIter = True
            print('results have not converged in the maximum ammount of iterations allowed')
    print ('Iterations have ended. Finished in ', n ,' itterations')
    return (MPhi, rho0, rhoS, E)



print('PG 0.2.7.C')

#Settings:
Test = False #set to True or False , Controls the plot functions

voltage = -80000 
density_constant = 1E-6 


roomHight = 1
roomWidth = 1

electrodePosX=0.5
electrodePosY=0.55
CondSize = 0.001 #in relation to whole

mesh_resolution = 100 

tol=1E-12

RoomTempC = 19
AirPressureC = 1
CondSurfCon = 1
RealCondSize = 1* CondSize

print('Generating: Mesh')

# Create mesh and define function space
space = Rectangle(Point(0,0), Point(roomWidth,roomHight))
cylinder = Circle(Point(electrodePosX, electrodePosY), CondSize)
domain = space - cylinder
mesh = generate_mesh(domain, mesh_resolution )

##########################################################################
# find the vertex closest to the center of the conductor

mesh_array=mesh.num_vertices()

coor = mesh.coordinates()   
disX = [(coor[index][0]-electrodePosX) for index in range(mesh_array)]

disY = [(coor[index][1] - electrodePosY) for index in range(mesh_array)]
dis = [np.sqrt(disX[index]**2+disY[index]**2) for index in range(len(disX))]

minDis = min(dis)
minDisIndex = dis.index(minDis)

print('Index of the Vertex with the min Distance: ', minDisIndex)
print('Distance to the Conductor center: ', dis[minDisIndex])
print('Coordinates of the Vertex with min Distance: ', coor[minDisIndex][0], coor[minDisIndex][1])

###########################################################################

if voltage <= 0:
    E0_Ec= 31*10**5
    K_Ec = 0.0308 

    #d_Ec = (273.15 + 20)/(273.15+RoomTempC)*(AirPressureC / 1013)
    d_Ec=1

    Ec=CondSurfCon * E0_Ec * d_Ec*(1+ K_Ec/(sqrt(d_Ec*RealCondSize)))

else :
    E0_Ec = 33.7 *10**5
    K_Ec = 0.024

    #d_Ec = (273.15 + 20)/(273.15+RoomTempC)*(AirPressureC / 1013)
    d_Ec=1

    Ec=CondSurfCon * E0_Ec * d_Ec*(1+ K_Ec/(sqrt(d_Ec*RealCondSize)))


V = VectorFunctionSpace(mesh, 'Lagrange', 2)  
Q = FunctionSpace(mesh, 'Lagrange', 2)

print('Calculating solution for Laplace Equation')

print ('using Neumann Border Conditions')


# Neumann Condition Setup
class Conductor(SubDomain):
    def inside(self, x, on_boundary):
        #return on_boundary and x[0]>0.48 and x[0]<0.51 and x[1]>0.54 and x[1]<0.56
        return (between(x[0], (electrodePosX-CondSize, electrodePosX+CondSize)) and between(x[1], (electrodePosY-CondSize, electrodePosY+CondSize)))


conductor = Conductor()

boundariesC = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundariesC.set_all(0)
conductor.mark(boundariesC, 1)

ds = ds(subdomain_data=boundariesC) 

#Ground Dirichlet Condition
ground = 'near(x[1], 0)'
bcp_ground = DirichletBC(Q, Constant(0), ground)

bcp = [bcp_ground]

#bcd = DirichletBC(Q, Constant(density_constant), conductor)
bcd = DirichletBC(Q, Constant(0), conductor)


print('Part: Laplace Equation = Start')
# Define variational problem
phi  = TrialFunction(Q)
v  = TestFunction(Q) 
f = Constant(0.0)


print("Ec = ", Ec)

# Laplace Equation for initial values
aL  = inner(grad(phi), grad(v))*dx
LL = f*v*dx - Ec * v *ds(1)

phi = Function(Q)

solve(aL == LL, phi, bcp) 

LaplacePhi=plot(phi, title='show laplace phi')
plt.colorbar(LaplacePhi)
plt.show()


##########################################################



if True:
    print('Laplace')
    phiL=phi

# control variables
LoopBreakPoint = 50
tol = 1E-10
minIter = 1
maxIter = 10

rho0 = None
rhoS = None
E = 0

RHO_CAHNGE_RATE = 1

while ((E-E0_Ec)/E0_Ec != tol):
    print("Calculate phi...")
    (phi, rho0, rhoS, E) = calculate_phi(phi, rhoS)
    print("Calculate phi, results:")
    print("phi = ", phi)
    print("rhoS = ", rhoS)
    print("Ec = ", E0_Ec)
    print("E = ", E)
    print("type of E = ", type(E))


    # Update rhoS
    if (lt(E,E0_Ec)):
        rhoS = rhoS - RHO_CAHNGE_RATE
    else:
        rhoS = rhoS + RHO_CAHNGE_RATE



file1 = File("phi-changing.pvd")
file1 << phi

file2 = File("rohS-changing.pvd")
file2 << rhoS

print(float(rho0))
#file3 = File("roh0-changing.pvd")
#file3 << rho0

#dolfin cant read pvd files
ChangingFile = HDF5File(MPI.comm_world,"static.h5","w")
ChangingFile.write(phi,"phi-changing")
ChangingFile.write(rhoS,"rhoS-changing")
#ChangingFile.write(rho0,"rho0-changing")
ChangingFile.close()

print('Result')
if True:
    plot(mesh, title='the mesh')
    plt.show()

    LaplacePhi=plot(phiL, title='show laplace phi')
    plt.colorbar(LaplacePhi)
    plt.show()

ResultE=plot(E, title='Result E')
plt.colorbar(ResultE) 
plt.show()

#rho=rho0+rhoS
ResultRho=plot(rhoS, title='Result rho')
plt.colorbar(ResultRho) 
plt.show()

ResultPhi=plot(phi, title='Result phi')
plt.colorbar(ResultPhi)
    
plt.show()


