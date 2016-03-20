
"""
Example from Brezzi's paper "A Family of Mimetic Finite 
Difference Methods on Polygonal and Polyhedral Mesh." 
"""

import sys

# Include relative path for mimpy library. 
sys.path.append("../../../")

import mimpy.mesh.voromesh as voromesh
import mimpy.mfd.mfd as mfd
import numpy as np
import random

res_mfd = mfd.MFD()
res_mfd.set_compute_diagonality(True)

method = 1
res_mfd.set_m_e_construction_method(method)

print "using method =>", method

import pickle

#Define the permeability function
def K(p):
    k = np.eye(3)
#    return k 
    x = p[0]
    y = p[1]
    z = p[2]
    k[0, 0] = 1.+y**2+z**2
    k[0, 1] = -x*y
    k[0, 2] = -x*z
    k[1, 0] = -x*y
    k[1, 1] = 1.+x**2+z**2
    k[1, 2] = -y*z
    k[2, 0] = -x*z
    k[2, 1] = -y*z
    k[2, 2] = 1.+x**2+y**2
    return k 
   
#The exact solution.
def u(p):
    x = p[0]
    y = p[1]
    z = p[2]
#    return x+2.*y-z
    return (x**3)*(y**2)*z+x*np.sin(2.*np.pi*x*y)*np.sin(2.*np.pi*y*z)*np.sin(2.*np.pi*z)

def dp_dx(p):
    x = p[0]
    y = p[1]
    z = p[2]
    return_value = 3.*x**2*y**2*z 
    return_value += 2.*np.pi*x*y*np.cos(2.*np.pi*x*y)*np.sin(2.*np.pi*z)*np.sin(2.*np.pi*y*z)
    return_value += np.sin(2.*np.pi*x*y)*np.sin(2.*np.pi*z)*np.sin(2.*np.pi*y*z)
    return return_value

def dp_dy(p):
    x = p[0]
    y = p[1]
    z = p[2]
    return_value = 2.*x**3*y*z+2.*np.pi*x*z*np.cos(2.*np.pi*y*z)*np.sin(2.*np.pi*x*y)*np.sin(2.*np.pi*z)  
    return_value += 2.*np.pi*x**2*np.cos(2.*np.pi*x*y)*np.sin(2.*np.pi*z)*np.sin(2.*np.pi*y*z)
    return return_value

def dp_dz(p):
    x = p[0]
    y = p[1]
    z = p[2]
    return_value = x**3*y**2+2.*np.pi*x*y*np.cos(2.*np.pi*y*z)*np.sin(2.*np.pi*x*y)*np.sin(2.*np.pi*z) 
    return_value += 2.*np.pi*x*np.cos(2.*np.pi*z)*np.sin(2.*np.pi*x*y)*np.sin(2.*np.pi*y*z)
    return return_value
    
def d2p_dx2(p):    
    x = p[0]
    y = p[1]
    z = p[2]
    return_value = 6.*x*y**2*z+4.*np.pi*y*np.cos(2.*np.pi*x*y)*np.sin(2.*np.pi*z)*np.sin(2.*np.pi*y*z)
    return_value -= 4.*np.pi**2*x*y**2*np.sin(2.*np.pi*x*y)*np.sin(2.*np.pi*z)*np.sin(2.*np.pi*y*z)
    return return_value
    
def d2p_dxy(p):
    x = p[0]
    y = p[1]
    z = p[2]
    return_value = 6.*x**2*y*z+4.*np.pi**2*x*y*z*np.cos(2.*np.pi*x*y)*np.cos(2.*np.pi*y*z)*np.sin(2.*np.pi*z)
    return_value += 2.*np.pi*z*np.cos(2.*np.pi*y*z)*np.sin(2.*np.pi*x*y)*np.sin(2.*np.pi*z)
    return_value += 4.*np.pi*x*np.cos(2.*np.pi*x*y)*np.sin(2.*np.pi*z)*np.sin(2.*np.pi*y*z)
    return_value -= 4.*np.pi**2*x**2*y*np.sin(2.*np.pi*x*y)*np.sin(2.*np.pi*z)*np.sin(2.*np.pi*y*z)
    return return_value

def d2p_dxz(p):
    x = p[0]
    y = p[1]
    z = p[2]
    return_value = 3.*x**2*y**2+4.*np.pi**2*x*y**2*np.cos(2.*np.pi*x*y)*np.cos(2.*np.pi*y*z)*np.sin(2.*np.pi*z)
    return_value += 2.*np.pi*y*np.cos(2.*np.pi*y*z)*np.sin(2.*np.pi*x*y)*np.sin(2.*np.pi*z)
    return_value += 4.*np.pi**2*x*y*np.cos(2.*np.pi*x*y)*np.cos(2.*np.pi*z)*np.sin(2.*np.pi*y*z)
    return_value += 2.*np.pi*np.cos(2.*np.pi*z)*np.sin(2.*np.pi*x*y)*np.sin(2.*np.pi*y*z)
    return return_value

def d2p_dy2(p):
    x = p[0]
    y = p[1]
    z = p[2]
    return_value = 2.* x**3*z+8.*np.pi**2*x**2*z*np.cos(2.*np.pi*x*y)*np.cos(2.*np.pi*y*z)*np.sin(2.*np.pi*z)
    return_value -= 4.*np.pi**2*x**3*np.sin(2.*np.pi*x*y)*np.sin(2.*np.pi*z)*np.sin(2.*np.pi*y*z)
    return_value -= 4.*np.pi**2*x*z**2*np.sin(2.*np.pi*x*y)*np.sin(2.*np.pi*z)*np.sin(2.*np.pi*y*z)
    return return_value

def d2p_dzy(p):
    x = p[0]
    y = p[1]
    z = p[2]
    return_value = 2*x**3*y+4.*np.pi**2*x*z*np.cos(2.*np.pi*z)*np.cos(2.*np.pi*y*z)*np.sin(2.*np.pi*x*y) 
    return_value += 4.*np.pi**2*x**2*y*np.cos(2.*np.pi*x*y)*np.cos(2.*np.pi*y*z)*np.sin(2.*np.pi*z) 
    return_value += 2.*np.pi*x*np.cos(2.*np.pi*y*z)*np.sin(2.*np.pi*x*y)*np.sin(2.*np.pi*z) 
    return_value += 4.*np.pi**2*x**2*np.cos(2.*np.pi*x*y)*np.cos(2.*np.pi*z)*np.sin(2.*np.pi*y*z) 
    return_value -= 4.*np.pi**2*x*y*z*np.sin(2.*np.pi*x*y)*np.sin(2.*np.pi*z)*np.sin(2.*np.pi*y*z)
    return return_value

def d2p_dz2(p):
    x = p[0]
    y = p[1]
    z = p[2]
    return_value = 8.*np.pi**2*x*y*np.cos(2.*np.pi*z)*np.cos(2.*np.pi*y*z)*np.sin(2.*np.pi*x*y)  
    return_value -= 4.*np.pi**2*x*np.sin(2.*np.pi*x*y)*np.sin(2.*np.pi*z)*np.sin(2.*np.pi*y*z)  
    return_value -= 4.*np.pi**2*x*y**2*np.sin(2.*np.pi*x*y)*np.sin(2.*np.pi*z)*np.sin(2.*np.pi*y*z)
    return return_value

def grad_u(p):
    x = p[0]
    y = p[1]
    z = p[2]
#    return -np.array([1., 2., -1.])
    return -np.array([(1.+y**2+z**2)*dp_dx(p)-x*y*dp_dy(p)-x*z*dp_dz(p), 
                      -x*y*dp_dx(p)+(1+x**2+z**2)*dp_dy(p)-y*z*dp_dz(p),
                      -x*z*dp_dx(p)-y*z*dp_dy(p)+(1+x**2+y**2)*dp_dz(p)])
                      
#The forcing function for the exact solution. 
def f(p):
    x = p[0]
    y = p[1]
    z = p[2]
#    return 0.
    return_value = (1.+y**2+z**2)*d2p_dx2(p)
    return_value += -y*dp_dy(p)-x*y*d2p_dxy(p)
    return_value += -z*dp_dz(p)-x*z*d2p_dxz(p)
    return_value += -x*dp_dx(p)-x*y*d2p_dxy(p)
    return_value += (1.+x**2+z**2)*d2p_dy2(p)
    return_value += -z*dp_dz(p)-y*z*d2p_dzy(p)
    return_value += -x*dp_dx(p)-x*z*d2p_dxz(p)
    return_value += -y*dp_dy(p)-y*z*d2p_dzy(p)
    return_value += (1.+x**2+y**2)*d2p_dz2(p)
    return -return_value

#The modification function is applied to the points of the mesh. 
#In this case no change is applied. 
def mod_function(p):
    return p
    return p +np.array([.02*np.sin(2.*np.pi*p[0])*np.sin(2.*np.pi*p[1])*np.sin(2.*np.pi*p[2]),
                         .02*np.sin(2.*np.pi*p[0])*np.sin(2.*np.pi*p[1])*np.sin(2.*np.pi*p[2]),
                         .02*np.sin(2.*np.pi*p[0])*np.sin(2.*np.pi*p[1])*np.sin(2.*np.pi*p[2])])

def make_random_mod_function(h):
    def mod(p, i, j, k):
        if .99>abs(p[0])>1.e-8 and \
                .99>abs(p[1])>1.e-8 and \
                .99>abs(p[2])>1.e-8:
            return p+np.array([random.random()*h*.2, 
                               random.random()*h*.0, 
                               random.random()*h*.0])
        return p
    return mod

def make_mod_function(h, ni, nj, nk):
    def mod(p, i, j, k):
        modification = np.array([0.,0.,0.])
        if j%2 == 0 and 0 < k < nk-1 :
            modification[2] += h*.25+h*random.random()*.2
        if i%2 == 0 and 0 < j < ni-1 :
            modification[1] += h*.0
        return p+modification
                          
    return mod


mesh_input = [(K, open("mesh4.vol")), 
              (K, open("mesh8.vol")), 
              (K, open("mesh16.vol")), 
              (K, open("mesh32.vol")), ]
#              (K, open("mesh64.vol")), ]



mesh_files = ['mesh4', 
              'mesh8', 
              'mesh16',
              'mesh32',]
#              'mesh64']

shift_parameter = .0

shifts = [1./4.,
          1./8.,
          1./16.,
          1./32.,]
#          1./64]
#          shift_parameter*1./32., ]
#          shift_parameter*1./64., ]
#          shift_parameter*1./128.,]
           
print "Centroid shifted by", shift_parameter

meshes = []
for (mesh_args, mesh_name) in zip(mesh_input, mesh_files):
    print mesh_name
    res_mesh = voromesh.VoroMesh()
    res_mesh.build_mesh(*mesh_args)   
    meshes.append(res_mesh)
#    cell_list = [0]*res_mesh.get_number_of_cells()
#    res_mesh.output_vtk_mesh(mesh_name+"blank", [map(lambda x:random.random(), cell_list)], ["random"])
#    pickle.dump(res_mesh, open(mesh_name, 'w'))

print "done building meshes" 


error_list = [[1./4.], 
              [1./8.], 
              [1./16.], 
              [1./32.],]
#              [1./64.], ]
#              [1./128.],]

      
for (index, (res_mesh, shift, mesh_name)) in enumerate(zip(meshes, shifts, mesh_files)):
    
    print "processing", mesh_name
#    res_mesh = pickle.load(open(mesh_name))
#    res_mesh.use_cell_shifted_centroid()
#    res_mesh.initialize_cell_shifted_centroid()

    cell_shifted = True
    point_shifted = True

    res_mesh.has_cell_shifted_centroid = cell_shifted
    res_mesh.has_face_shifted_centroid = point_shifted
    
    print "cell shifted = ", cell_shifted
    print "piont shfited = ", point_shifted

    res_mesh.apply_dirichlet_from_function(0, lambda p:u(p))
    res_mesh.apply_dirichlet_from_function(1, lambda p:u(p))
    res_mesh.apply_dirichlet_from_function(2, lambda p:u(p))
    res_mesh.apply_dirichlet_from_function(3, lambda p:u(p))
    res_mesh.apply_dirichlet_from_function(4, lambda p:u(p))
    res_mesh.apply_dirichlet_from_function(5, lambda p:u(p))

    #Apply the forcing function f. 
    res_mesh.apply_forcing_from_function(f)
    
    #res_mesh.apply_forcing_from_grad(grad_u, f)
    
    #Connect the MFD instance to the new mesh. 
    res_mfd.set_mesh(res_mesh)

    res_mfd.check_m_e = True

    print "build lhs"
    #Build the LHS and RHS. 
    #res_mfd.build_lhs_divided()
    res_mfd.build_lhs_petsc()
    #res_mfd.build_lhs()
    print "build rhs"
    res_mfd.build_rhs()
    
    #pickle.dump(res_mfd.lhs, open(mesh_name+"_COO", 'w'))
    
    #Solve the linear system. 
    #res_mfd.solve(solver = 'gmres')
    #res_mfd.solve_divided()
    res_mfd.solve_petsc()

    #res_mesh.output_vtk_mesh(mesh_name, 
    #                         [res_mfd.get_pressure_solution(), 
    #                          res_mfd.get_analytical_pressure_solution(u), 
    #                          res_mfd.get_diagonality(), ],
    #                         ["MFDPressure", "AnalyticalPressure", "ORTHO"])

    error_list[index] += [res_mfd.compute_l2_error_pressure(u),
                          res_mfd.compute_l2_error_pressure(u, quadrature_method = 1),
                          res_mfd.compute_l2_error_velocity(grad_u), 
                          res_mfd.compute_l2_error_velocity(grad_u, 1),]

print error_list[0][0], "&" , "%.4e"%error_list[0][1], "&" ,
print "---", "&" ,"%.4e"%error_list[0][2],"&" , "---", "&","%.4e"%error_list[0][3], "&" ,"---", 
print "%.4e"%error_list[0][4], "---"
for i in range(1, len(mesh_files)):
    print "%.4e &"% error_list[i][0],
    rate_1 = ((np.log(error_list[i][1])-np.log(error_list[i-1][1]))/
              (np.log(error_list[i][0])-np.log(error_list[i-1][0])))
    print "%.4e & %.4f &" % (error_list[i][1], rate_1),  
    rate_2 = ((np.log(error_list[i][2])-np.log(error_list[i-1][2]))/
              (np.log(error_list[i][0])-np.log(error_list[i-1][0])))
    print "%.4e & %.4f &" % (error_list[i][2],rate_2), 
    rate_3 = ((np.log(error_list[i][3])-np.log(error_list[i-1][3]))/
              (np.log(error_list[i][0])-np.log(error_list[i-1][0])))
    print "%.4e & %.4f" % (error_list[i][3],rate_3),
    rate_4 = ((np.log(error_list[i][4])-np.log(error_list[i-1][4]))/
              (np.log(error_list[i][0])-np.log(error_list[i-1][0])))
    print "%.4e & %.4f" % (error_list[i][4],rate_4)
