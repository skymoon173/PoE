# -*- coding: utf-8 -*-
# coding=utf-8
# coding: utf-8

"""
Rx(θ) = [[1    0            0      0] 
             [0 cosθ -sinθ   0]
			 [0 sinθ cosθ 0]
			 [0    0            0      1]]
			 
Ry(θ) = [[cosθ   0     sinθ      0] 
             [0             1        0             0]
			 [-sinθ  0     cosθ      0]
			 [0             0        0             1]]

Rz(θ) = [[cosθ   -sinθ    0     0] 
    	 [sinθ    cosθ    0     0]
		 [0        0      1     0]
		 [0        0      0     1]]


Tn = [[cos(θ_n)                         -sin(θ_n)              0                  a_n-1      ]
	  [sin(θ_n)*cos(α_n-1)   cos(θ_n)*cos(α_n-1) -sin(α_n-1)  -d_n*sin(α_n-1)]
	  [sin(θ_n)*sin(α_n-1)   cos(θ_n)*sin(α_n-1)  cos(α_n-1)   d_n*cos(α_n-1)]
	  [    0                                       0                   0                    1        ]]

			 
T06 = [[normal_x    orientation_x   approach_x    Px]
	   [normal_y    orientation_y   approach_y    Py]
	   [normal_z    orientation_z   approach_z    Pz]
	   [   0      		  0              0         1]]

Euler(ϕ, θ, ψ) = Rz(ϕ)*Ry(θ)*Rz(ψ)
RPY(ϕ, θ, ψ) = Rz(ϕ)*Ry(θ)*Rx(ψ)
                     =[[cosϕ*cosθ  -sinϕ*cosψ+cosϕ*cosθ*sinψ     sinϕ*sinψ+cosϕ*sinθ*cosψ]
	 				   [sinϕ*cosθ   cosϕ*cosψ+sinϕ*sinθ*sinψ    -cosϕ*sinψ+sinϕ*sinθ*cosψ]
	  				   [-sinθ               cosθ*sinψ                   cosθ*cosψ        ]]


"""
import sympy
import numpy as np

"""RPY角旋转矩阵公式"""
cosϕ = sympy.Symbol('cosϕ')
sinϕ = sympy.Symbol('sinϕ') 
cosθ = sympy.Symbol('cosθ')
sinθ = sympy.Symbol('sinθ') 
cosψ = sympy.Symbol('cosψ')
sinψ = sympy.Symbol('sinψ') 

Rz = sympy.Matrix([[cosϕ, -sinϕ,  0, 0], 
            	   [sinϕ,  cosϕ,  0 ,0],
				   [0,       0 ,  1, 0],
				   [0,       0,   0 ,1]])

Ry = sympy.Matrix([[ cosθ,   0,     sinθ,    0], 
            	   [  0,     1,      0,      0],
			       [-sinθ,   0,     cosθ,    0],
			 	   [  0,     0,      0,      1]])

Rx = sympy.Matrix([[1,  0,     0,   0], 
             	   [0 ,cosψ, -sinψ, 0],
				   [0, sinψ,  cosψ, 0],
			 	   [0,  0,     0,   1]])

'''RPY角运算公式'''
RPY_R = Rz*Ry*Rx

RPY_Rotation =  sympy.Matrix([[cosθ*cosϕ, -cosψ*sinϕ + cosϕ*sinθ*sinψ, cosψ*cosϕ*sinθ + sinψ*sinϕ, 0], 
						 	 [cosθ*sinϕ, cosψ*cosϕ + sinθ*sinψ*sinϕ, cosψ*sinθ*sinϕ - cosϕ*sinψ, 0], 
						     [-sinθ, cosθ*sinψ, cosθ*cosψ, 0], 
						 	 [0, 0, 0, 1]])


"""DH坐标变换公式"""
θ1 = sympy.Symbol('θ1')
θ2 = sympy.Symbol('θ2') 
θ3 = 0
θ4 = sympy.Symbol('θ4') 

cosθ1 = sympy.Symbol('cosθ1')
cosθ2 = sympy.Symbol('cosθ2') 
cosθ3 = 1
cosθ4 = sympy.Symbol('cosθ4') 

sinθ1 = sympy.Symbol('sinθ1') 
sinθ2 = sympy.Symbol('sinθ2') 
sinθ3 = 0
sinθ4 = sympy.Symbol('sinθ4') 

sinα0 = 0
sinα1 = 0
sinα2 = 0
sinα3 = 0

cosα0 = 1
cosα1 = 1
cosα2 = 1
cosα3 = 1

a0  = 0
a1  = 325
a2  = 275
a3  = 0

J3 = sympy.Symbol('J3') 

d1  = 205
d2  = 0
d3  = J3-55
d4  = 0

T01 = sympy.Matrix([[cosθ1,-sinθ1,0,a0],[sinθ1*cosα0,cosθ1*cosα0,-sinα0,-d1*sinα0],[sinθ1*sinα0,cosθ1*sinα0,cosα0,d1*cosα0], [0,0,0,1]])
T12 = sympy.Matrix([[cosθ2,-sinθ2,0,a1],[sinθ2*cosα1,cosθ2*cosα1,-sinα1,-d2*sinα1],[sinθ2*sinα1,cosθ2*sinα1,cosα1,d2*cosα1], [0,0,0,1]])
T23 = sympy.Matrix([[cosθ3,-sinθ3,0,a2],[sinθ3*cosα2,cosθ3*cosα2,-sinα2,-d3*sinα2],[sinθ3*sinα2,cosθ3*sinα2,cosα2,d3*cosα2], [0,0,0,1]])
T34 = sympy.Matrix([[cosθ4,-sinθ4,0,a3],[sinθ4*cosα3,cosθ4*cosα3,-sinα3,-d4*sinα3],[sinθ4*sinα3,cosθ4*sinα3,cosα3,d4*cosα3], [0,0,0,1]])





'''
T02 = T01*T12
T03 = T01*T12*T23

T05 = T01*T12*T23*T34*T45
T06 = T01*T12*T23*T34*T45*T56	
  
T02_Inverse = T02**(-1)
'''

ax = sympy.Symbol('ax') 
ay = sympy.Symbol('ay')
az = sympy.Symbol('az') 
ox = sympy.Symbol('ox') 
oy = sympy.Symbol('oy')
oz = sympy.Symbol('oz') 

nx = sympy.Symbol('nx') 
ny = sympy.Symbol('ny')
nz = sympy.Symbol('nz') 

Px = sympy.Symbol('Px') 
Py = sympy.Symbol('Py')
Pz = sympy.Symbol('Pz') 

'''计算
T06 = T03 * T14
T03_Inverse * T06 = T14
T14_Inverse * T03_Inverse * T06 = E
'''
T01_Inverse =  sympy.Matrix([[cosθ1/(cosθ1**2 + sinθ1**2), sinθ1/(cosθ1**2 + sinθ1**2), 0, 0], [-sinθ1/(cosθ1**2 + sinθ1**2), cosθ1/(cosθ1**2 + sinθ1**2), 0, 0], [0, 0, 1, -205], [0, 0, 0, 1]])
'''
T04 = sympy.Matrix([[cosθ*cosϕ, -cosψ*sinϕ + cosϕ*sinθ*sinψ, cosψ*cosϕ*sinθ + sinψ*sinϕ, Px],[cosθ*sinϕ, cosψ*cosϕ + sinθ*sinψ*sinϕ, cosψ*sinθ*sinϕ - cosϕ*sinψ, Py], [-sinθ, cosθ*sinψ, cosθ*cosψ, Pz],[0,  0,  0,  1 ]])
T01_Inverse_T04 =  sympy.Matrix([[cosθ*cosθ1*cosϕ/(cosθ1**2 + sinθ1**2) + cosθ*sinθ1*sinϕ/(cosθ1**2 + sinθ1**2), cosθ1*(-cosψ*sinϕ + cosϕ*sinθ*sinψ)/(cosθ1**2 + sinθ1**2) + sinθ1*(cosψ*cosϕ + sinθ*sinψ*sinϕ)/(cosθ1**2 + sinθ1**2), cosθ1*(cosψ*cosϕ*sinθ + sinψ*sinϕ)/(cosθ1**2 + sinθ1**2) + sinθ1*(cosψ*sinθ*sinϕ - cosϕ*sinψ)/(cosθ1**2 + sinθ1**2), Px*cosθ1/(cosθ1**2 + sinθ1**2) + Py*sinθ1/(cosθ1**2 + sinθ1**2)], [cosθ*cosθ1*sinϕ/(cosθ1**2 + sinθ1**2) - cosθ*cosϕ*sinθ1/(cosθ1**2 + sinθ1**2), cosθ1*(cosψ*cosϕ + sinθ*sinψ*sinϕ)/(cosθ1**2 + sinθ1**2) - sinθ1*(-cosψ*sinϕ + cosϕ*sinθ*sinψ)/(cosθ1**2 + sinθ1**2), cosθ1*(cosψ*sinθ*sinϕ - cosϕ*sinψ)/(cosθ1**2 + sinθ1**2) - sinθ1*(cosψ*cosϕ*sinθ + sinψ*sinϕ)/(cosθ1**2 + sinθ1**2), -Px*sinθ1/(cosθ1**2 + sinθ1**2) + Py*cosθ1/(cosθ1**2 + sinθ1**2)], [-sinθ, cosθ*sinψ, cosθ*cosψ, Pz - 205], [0, 0, 0, 1]])
T01_Inverse_simplize_T04 = sympy.Matrix([[cosθ*cosθ1*cosϕ + cosθ*sinθ1*sinϕ, cosθ1*(-cosψ*sinϕ + cosϕ*sinθ*sinψ) + sinθ1*(cosψ*cosϕ + sinθ*sinψ*sinϕ), cosθ1*(cosψ*cosϕ*sinθ + sinψ*sinϕ) + sinθ1*(cosψ*sinθ*sinϕ - cosϕ*sinψ), Px*cosθ1 + Py*sinθ1], 
[cosθ*cosθ1*sinϕ - cosθ*cosϕ*sinθ1, cosθ1*(cosψ*cosϕ + sinθ*sinψ*sinϕ) - sinθ1*(-cosψ*sinϕ + cosϕ*sinθ*sinψ), cosθ1*(cosψ*sinθ*sinϕ - cosϕ*sinψ) - sinθ1*(cosψ*cosϕ*sinθ + sinψ*sinϕ), -Px*sinθ1 + Py*cosθ1], 
[-sinθ, cosθ*sinψ, cosθ*cosψ, Pz - 205], 
[0, 0, 0, 1]])

'''


'''
expr2 = [T01_Inverse_simplize_T04[0,0]-T14[0,0], T01_Inverse_simplize_T04[0,1]-T14[0,1],T01_Inverse_simplize_T04[0,2]-T14[0,2],T01_Inverse_simplize_T04[0,3]-T14[0,3],T01_Inverse_simplize_T04[1,0]-T14[1,0],T01_Inverse_simplize_T04[1,1]-T14[1,1],T01_Inverse_simplize_T04[1,2]-T14[1,2],T01_Inverse_simplize_T04[1,3]-T14[1,3],T01_Inverse_simplize_T04[2,0]-T14[2,0],T01_Inverse_simplize_T04[2,1]-T14[2,1],T01_Inverse_simplize_T04[2,2]-T14[2,2],T01_Inverse_simplize_T04[2,3]-T14[2,3]]
r2 = sympy.solve(expr2, [cosθ1,cosθ2,cosθ4,sinθ1,sinθ2,sinθ4])
'''




ω = sympy.Symbol('ω')
ω1 = sympy.Symbol('ω1')
ω2 = sympy.Symbol('ω2')
ω3 = sympy.Symbol('ω3')

v = sympy.Symbol('v')
v1 = sympy.Symbol('v1')
v2 = sympy.Symbol('v2')
v3 = sympy.Symbol('v3')

θ = sympy.Symbol('θ')


ω  = sympy.Matrix([[0, -ω3, ω2], 
             	   [ω3, 0, -ω1],
			 	   [-ω2, ω1, 0]])



'''Rodrigues’s formula
Rot = I + sinθ*ω + (1-cosθ)*ω*ω
'''


class PoE:

	def S(ω1,ω2,ω3,v1,v2,v3):
		S = np.random.random(6)
		S = S.reshape(6, 1)

		S[0,0] = ω1              
		S[1,0] = ω2              
		S[2,0] = ω3  
		S[3,0] = v1              
		S[4,0] = v2              
		S[5,0] = v3 
		return 

	def e_ωθ(ω1,ω2,ω3,θ):
		I = np.random.random(9)
		I = I.reshape(3, 3)
		I[0,0] = 1               
		I[0,1] = 0                   
		I[0,2] = 0         
		I[1,0] = 0               
		I[1,1] = 1                   
		I[1,2] = 0    
		I[2,0] = 0               
		I[2,1] = 0                   
		I[2,2] = 1            
		
		ω = np.random.random(9)
		ω = ω.reshape(3, 3)

		ω[0,0] = 0               
		ω[0,1] = -ω3                   
		ω[0,2] = ω2        
		ω[1,0] = ω3              
		ω[1,1] = 0                  
		ω[1,2] = -ω1    
		ω[2,0] = -ω2               
		ω[2,1] = ω1                   
		ω[2,2] = 0  

		ω_2 = np.dot(ω,ω)
		R = I + np.dot(np.sin(θ),ω) + np.dot((1-np.cos(θ)),ω_2)

		return(R)

	def e_Bθ(ω1,ω2,ω3,v1,v2,v3,θ):
		I = np.random.random(9)
		I = I.reshape(3, 3)
		I[0,0] = 1               
		I[0,1] = 0                   
		I[0,2] = 0         
		I[1,0] = 0               
		I[1,1] = 1                   
		I[1,2] = 0    
		I[2,0] = 0               
		I[2,1] = 0                   
		I[2,2] = 1   

		ω = np.random.random(9)
		ω = ω.reshape(3, 3)

		ω[0,0] = 0               
		ω[0,1] = -ω3                   
		ω[0,2] = ω2        
		ω[1,0] = ω3              
		ω[1,1] = 0                  
		ω[1,2] = -ω1    
		ω[2,0] = -ω2               
		ω[2,1] = ω1                   
		ω[2,2] = 0  

		v = np.random.random(3)
		v = v.reshape(3, 1)
		v[0,0] = v1                     
		v[1,0] = v2                
		v[2,0] = v3              

		e_Bθ = np.random.random(16)
		e_Bθ = e_Bθ.reshape(4, 4)

		R = PoE.e_ωθ(ω1,ω2,ω3,θ)

		e_Bθ[0,0] = R[0,0]
		e_Bθ[1,0] = R[1,0]
		e_Bθ[2,0] = R[2,0]
		e_Bθ[3,0] = 0

		e_Bθ[0,1] = R[0,1]
		e_Bθ[1,1] = R[1,1]
		e_Bθ[2,1] = R[2,1]
		e_Bθ[3,1] = 0

		e_Bθ[0,2] = R[0,2]
		e_Bθ[1,2] = R[1,2]
		e_Bθ[2,2] = R[2,2]
		e_Bθ[3,2] = 0

		e_Bθ = e_Bθ.reshape(4, 4)
		V = np.dot(np.dot(I,θ)+np.dot((1-np.cos(θ)),ω)+np.dot((θ-np.sin(θ)),np.dot(ω,ω)),v)

		e_Bθ[0,3] = V[0]
		e_Bθ[1,3] = V[1]
		e_Bθ[2,3] = V[2]
		e_Bθ[3,3] = 1
		
		return e_Bθ

	def AdT(ω1,ω2,ω3,v1,v2,v3,θ):
		T = PoE.e_Bθ(ω1,ω2,ω3,v1,v2,v3,θ)

		R = np.random.random(9)
		R = R.reshape(3, 3)

		R[0,0] = T[0,0]
		R[1,0] = T[1,0]
		R[2,0] = T[2,0]

		R[0,1] = T[0,1]
		R[1,1] = T[1,1]
		R[2,1] = T[2,1]

		R[0,2] = T[0,2]
		R[1,2] = T[1,2]
		R[2,2] = T[2,2]

		p = np.random.random(3)
		p = p.reshape(3, 1)

		p[0,0] = T[0,3]
		p[1,0] = T[1,3]
		p[2,0] = T[2,3]

		P = np.random.random(9)
		P = P.reshape(3, 3)

		P[0,0] = 0               
		P[0,1] = -p[2,0]                  
		P[0,2] = p[1,0]        
		P[1,0] = p[2,0]             
		P[1,1] = 0                  
		P[1,2] = -p[0,0]   
		P[2,0] = -p[1,0]              
		P[2,1] = p[0,0]                  
		P[2,2] = 0  
		
		PR = np.dot(P,R)

		AdT = np.random.random(36)
		AdT = AdT.reshape(6, 6)

		AdT[0,0] = R[0,0]
		AdT[1,0] = R[1,0]
		AdT[2,0] = R[2,0]
		AdT[0,1] = R[0,1]
		AdT[1,1] = R[1,1]
		AdT[2,1] = R[2,1]
		AdT[0,2] = R[0,2]
		AdT[1,2] = R[1,2]
		AdT[2,2] = R[2,2]

		AdT[3,0] = PR[0,0]
		AdT[4,0] = PR[1,0]
		AdT[5,0] = PR[2,0]
		AdT[3,1] = PR[0,1]
		AdT[4,1] = PR[1,1]
		AdT[5,1] = PR[2,1]
		AdT[3,2] = PR[0,2]
		AdT[4,2] = PR[1,2]
		AdT[5,2] = PR[2,2]

		AdT[0,3] = 0
		AdT[1,3] = 0
		AdT[2,3] = 0
		AdT[0,4] = 0
		AdT[1,4] = 0
		AdT[2,4] = 0
		AdT[0,5] = 0
		AdT[1,5] = 0
		AdT[2,5] = 0

		AdT[3,3] = R[0,0]
		AdT[4,3] = R[1,0]
		AdT[5,3] = R[2,0]
		AdT[3,4] = R[0,1]
		AdT[4,4] = R[1,1]
		AdT[5,4] = R[2,1]
		AdT[3,5] = R[0,2]
		AdT[4,5] = R[1,2]
		AdT[5,5] = R[2,2]

		return AdT


def MZ07_Forward_Kinematics_PoEModel(θ1,θ2,θ3,θ4,θ5,θ6):

	M = np.array([[1,0,0,2],
				  [0,1,0,0],
				  [0,0,1,0],
				  [0,0,0,1]])


	S1 = np.array([[0],[0],[1],[0],[2],[0]])
	S2 = np.array([[0],[0],[1],[0],[1],[0]])

	S1_T = S1.T[0]
	S2_T = S2.T[0]


	e_Bθ1 = PoE.e_Bθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ1*np.pi/180)
	e_Bθ2 = PoE.e_Bθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ2*np.pi/180)


	T12 = np.dot(e_Bθ1, e_Bθ2)
	T = np.dot(M,T12) 


	AdT1 = PoE.AdT(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ1*np.pi/180)
	AdT2 = PoE.AdT(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ2*np.pi/180)

	AdT12 = np.dot(AdT1,AdT2)

	Js1 = S1
	Js2 = np.dot(AdT1,S2)

	Js = np.random.random(12)
	Js = Js.reshape(6, 2)

	Js[0,0] = Js1[0][0]
	Js[1,0] = Js1[1][0]
	Js[2,0] = Js1[2][0]
	Js[3,0] = Js1[3][0]
	Js[4,0] = Js1[4][0]
	Js[5,0] = Js1[5][0]

	Js[0,1] = Js2[0][0]
	Js[1,1] = Js2[1][0]
	Js[2,1] = Js2[2][0]
	Js[3,1] = Js2[3][0]
	Js[4,1] = Js2[4][0]
	Js[5,1] = Js2[5][0]


	'''求矩阵的伪逆(瘦高型)'''	
	Js_T = Js.T	
	Js_pseudoInverse = np.dot(np.linalg.inv(np.dot(Js_T,Js)),Js_T)


	print('Js_pseudoInverse =', Js_pseudoInverse)

	Euler_phi = np.arctan(ay/ax)
	Euler_theta = np.arctan((np.cos(Euler_phi)*ax+np.sin(Euler_phi)*ay)/az)
	#Euler_theta = np.arctan((ax+ay)/((np.cos(Euler_phi)+np.sin(Euler_phi))*az))
	#Euler_theta = np.arctan(ay/(az*np.sin(Euler_phi)))
	Euler_psi = np.arctan((-np.sin(Euler_phi)*nx+np.cos(Euler_phi)*ny)/(-np.sin(Euler_phi)*ox+np.cos(Euler_phi)*oy))
	#Euler_psi = np.arctan(-oz/nz)
	
	RPY_phi = np.arctan(ny/nx)
	#RPY_theta = np.arctan(-nz*(np.cos(RPY_phi)+np.sin(RPY_phi))/(nx+ny))
	RPY_theta = np.arctan(-nz/(np.cos(RPY_phi)*nx+np.sin(RPY_phi)*ny))
	#RPY_theta = np.arctan(-np.sin(RPY_phi)*nz/ny)*180/np.pi
	#RPY_theta = np.arcsin(-nz)
	RPY_psi = np.arctan((np.sin(RPY_phi)*ax-np.cos(RPY_phi)*ay)/(-np.sin(RPY_phi)*ox+np.cos(RPY_phi)*oy))
	#RPY_psi = np.arctan(oz/az)
	print('PoE = ', T)
	return [px, py, pz, RPY_phi*180/np.pi, RPY_theta*180/np.pi, RPY_psi*180/np.pi, Euler_phi*180/np.pi,Euler_theta*180/np.pi,Euler_psi*180/np.pi]


'''
e_Bθ1 = PoE.e_Sθ(0,0,-1,0,-50-330-450,0,θ1)
e_Bθ2 = PoE.e_Sθ(0,1,0,,0,0,θ2)
e_Bθ3 = PoE.e_Sθ(0,0,-1,-19,0,0,θ3)
e_Bθ4 = PoE.e_Sθ(0,0,0,0,0,-1,θ4)


T12_2 = np.dot(e_Bθ1, e_Bθ2)
T23_2 = np.dot(T12_2, e_Bθ3)
T3_2 = np.dot(M, T23_2)
print(T3_2)
'''

def EC06_Forward_Kinematics_PoEModel(θ1,θ2,θ3,θ4):

	M = np.random.random(16)
	M = M.reshape(4, 4)
	M[0,0] = 1
	M[1,0] = 0
	M[2,0] = 0
	M[3,0] = 0

	M[0,1] = 0
	M[1,1] = 1
	M[2,1] = 0
	M[3,1] = 0

	M[0,2] = 0
	M[1,2] = 0
	M[2,2] = 1
	M[3,2] = 0

	M[0,3] = 325+275
	M[1,3] = 0
	M[2,3] = 150
	M[3,3] = 1

	e_Sθ1 = PoE.e_Sθ(0,0,1,0,0,0,θ1*np.pi/180)
	e_Sθ2 = PoE.e_Sθ(0,0,1,0,-325,0,θ2*np.pi/180)
	e_Sθ3 = PoE.e_Sθ(0,0,0,0,0,1,θ3)
	e_Sθ4 = PoE.e_Sθ(0,0,1,0,-325-275,0,θ4*np.pi/180)



	T12 = np.dot(e_Sθ1, e_Sθ2)
	T13 = np.dot(T12, e_Sθ3)
	T14 = np.dot(T13, e_Sθ4)
	
	

	T = np.dot(T14, M) 

	nx = T[0,0]
	ny = T[1,0]
	nz = T[2,0]
	ox = T[0,1]
	oy = T[1,1]
	oz = T[2,1]
	ax = T[0,2]
	ay = T[1,2]
	az = T[2,2]
	px = T[0,3]
	py = T[1,3]
	pz = T[2,3]

	Euler_phi = np.arctan(ay/ax)
	Euler_theta = np.arctan((np.cos(Euler_phi)*ax+np.sin(Euler_phi)*ay)/az)
	#Euler_theta = np.arctan((ax+ay)/((np.cos(Euler_phi)+np.sin(Euler_phi))*az))
	#Euler_theta = np.arctan(ay/(az*np.sin(Euler_phi)))
	Euler_psi = np.arctan((-np.sin(Euler_phi)*nx+np.cos(Euler_phi)*ny)/(-np.sin(Euler_phi)*ox+np.cos(Euler_phi)*oy))
	#Euler_psi = np.arctan(-oz/nz)
	
	RPY_phi = np.arctan(ny/nx)
	#RPY_theta = np.arctan(-nz*(np.cos(RPY_phi)+np.sin(RPY_phi))/(nx+ny))
	RPY_theta = np.arctan(-nz/(np.cos(RPY_phi)*nx+np.sin(RPY_phi)*ny))
	#RPY_theta = np.arctan(-np.sin(RPY_phi)*nz/ny)*180/np.pi
	#RPY_theta = np.arcsin(-nz)
	RPY_psi = np.arctan((np.sin(RPY_phi)*ax-np.cos(RPY_phi)*ay)/(-np.sin(RPY_phi)*ox+np.cos(RPY_phi)*oy))
	#RPY_psi = np.arctan(oz/az)
	return [px, py, pz, RPY_phi*180/np.pi, RPY_theta*180/np.pi, RPY_psi*180/np.pi, Euler_phi*180/np.pi,Euler_theta*180/np.pi,Euler_psi*180/np.pi]


# -*- coding: utf-8 -*-
# coding=utf-8
# coding: utf-8

"""
Rx(θ) = [[1    0            0      0] 
             [0 cosθ -sinθ   0]
			 [0 sinθ cosθ 0]
			 [0    0            0      1]]
			 
Ry(θ) = [[cosθ   0     sinθ      0] 
             [0             1        0             0]
			 [-sinθ  0     cosθ      0]
			 [0             0        0             1]]

Rz(θ) = [[cosθ   -sinθ    0     0] 
    	 [sinθ    cosθ    0     0]
		 [0        0      1     0]
		 [0        0      0     1]]


Tn = [[cos(θ_n)                         -sin(θ_n)              0                  a_n-1      ]
	  [sin(θ_n)*cos(α_n-1)   cos(θ_n)*cos(α_n-1) -sin(α_n-1)  -d_n*sin(α_n-1)]
	  [sin(θ_n)*sin(α_n-1)   cos(θ_n)*sin(α_n-1)  cos(α_n-1)   d_n*cos(α_n-1)]
	  [    0                                       0                   0                    1        ]]

			 
T06 = [[normal_x    orientation_x   approach_x    Px]
	   [normal_y    orientation_y   approach_y    Py]
	   [normal_z    orientation_z   approach_z    Pz]
	   [   0      		  0              0         1]]

Euler(ϕ, θ, ψ) = Rz(ϕ)*Ry(θ)*Rz(ψ)
RPY(ϕ, θ, ψ) = Rz(ϕ)*Ry(θ)*Rx(ψ)
                     =[[cosϕ*cosθ  -sinϕ*cosψ+cosϕ*cosθ*sinψ     sinϕ*sinψ+cosϕ*cosθ*cosψ]
	 				   [sinϕ*cosθ   cosϕ*cosψ+sinϕ*sinθ*sinψ    -cosϕ*sinψ+sinϕ*sinθ*cosψ]
	  				   [-sinθ               cosθ*sinψ                   cosθ*cosψ        ]]


"""
import sympy
import numpy as np

"""RPY角旋转矩阵公式"""
cosϕ = sympy.Symbol('cosϕ')
sinϕ = sympy.Symbol('sinϕ') 
cosθ = sympy.Symbol('cosθ')
sinθ = sympy.Symbol('sinθ') 
cosψ = sympy.Symbol('cosψ')
sinψ = sympy.Symbol('sinψ') 

Rz = sympy.Matrix([[cosϕ, -sinϕ,  0, 0], 
            	   [sinϕ,  cosϕ,  0 ,0],
				   [0,       0 ,  1, 0],
				   [0,       0,   0 ,1]])

Ry = sympy.Matrix([[ cosθ,   0,     sinθ,    0], 
            	   [  0,     1,      0,      0],
			       [-sinθ,   0,     cosθ,    0],
			 	   [  0,     0,      0,      1]])

Rx = sympy.Matrix([[1,  0,     0,   0], 
             	   [0 ,cosψ, -sinψ, 0],
				   [0, sinψ,  cosψ, 0],
			 	   [0,  0,     0,   1]])

'''RPY角运算公式'''
RPY_R = Rz*Ry*Rx

RPY_Rotation =  sympy.Matrix([[cosθ*cosϕ, -cosψ*sinϕ + cosϕ*sinθ*sinψ, cosψ*cosϕ*sinθ + sinψ*sinϕ, 0], 
						 	 [cosθ*sinϕ, cosψ*cosϕ + sinθ*sinψ*sinϕ, cosψ*sinθ*sinϕ - cosϕ*sinψ, 0], 
						     [-sinθ, cosθ*sinψ, cosθ*cosψ, 0], 
						 	 [0, 0, 0, 1]])


"""DH坐标变换公式"""
θ1 = sympy.Symbol('θ1')
θ2 = sympy.Symbol('θ2') 
θ3 = 0
θ4 = sympy.Symbol('θ4')
θ1θ2θ4 = sympy.Symbol('θ1θ2θ4') 

cosθ1 = sympy.Symbol('cosθ1')
cosθ2 = sympy.Symbol('cosθ2') 
cosθ3 = 1
cosθ4 = sympy.Symbol('cosθ4') 
cosθ1θ2 = sympy.Symbol('cosθ1θ2') 
cosθ1θ2θ4 = sympy.Symbol('cosθ1θ2θ4') 

sinθ1 = sympy.Symbol('sinθ1') 
sinθ2 = sympy.Symbol('sinθ2') 
sinθ3 = 0
sinθ4 = sympy.Symbol('sinθ4') 
sinθ1θ2 = sympy.Symbol('sinθ1θ2') 
sinθ1θ2θ4 = sympy.Symbol('sinθ1θ2θ4') 


sinα0 = 0
sinα1 = 0
sinα2 = 0
sinα3 = 0

cosα0 = 1
cosα1 = 1
cosα2 = 1
cosα3 = 1

a0  = 0
a1  = 325
a2  = 275
a3  = 0

J3 = sympy.Symbol('J3') 

v1 = sympy.Symbol('v1') 
v2 = sympy.Symbol('v2') 
v3 = sympy.Symbol('v3') 

px = sympy.Symbol('px') 
py = sympy.Symbol('py') 
pz = sympy.Symbol('pz') 

Px = sympy.Symbol('Px') 
Py = sympy.Symbol('Py') 
Pz = sympy.Symbol('Pz') 

ox = sympy.Symbol('ox') 
oy = sympy.Symbol('oy') 
V1 = sympy.Symbol('V1') 
V2 = sympy.Symbol('V2') 


d1  = 205
d2  = 0
d3  = J3-55
d4  = 0

T01 = sympy.Matrix([[cosθ1,-sinθ1,0,a0],[sinθ1*cosα0,cosθ1*cosα0,-sinα0,-d1*sinα0],[sinθ1*sinα0,cosθ1*sinα0,cosα0,d1*cosα0], [0,0,0,1]])
T12 = sympy.Matrix([[cosθ2,-sinθ2,0,a1],[sinθ2*cosα1,cosθ2*cosα1,-sinα1,-d2*sinα1],[sinθ2*sinα1,cosθ2*sinα1,cosα1,d2*cosα1], [0,0,0,1]])
T23 = sympy.Matrix([[cosθ3,-sinθ3,0,a2],[sinθ3*cosα2,cosθ3*cosα2,-sinα2,-d3*sinα2],[sinθ3*sinα2,cosθ3*sinα2,cosα2,d3*cosα2], [0,0,0,1]])
T34 = sympy.Matrix([[cosθ4,-sinθ4,0,a3],[sinθ4*cosα3,cosθ4*cosα3,-sinα3,-d4*sinα3],[sinθ4*sinα3,cosθ4*sinα3,cosα3,d4*cosα3], [0,0,0,1]])

T02 = np.dot(T01,T12)
T03 = np.dot(T02,T23)
T04 = np.dot(T03,T34)

'''
T04 =  [[cosθ4*(cosθ1*cosθ2 - sinθ1*sinθ2) + sinθ4*(-cosθ1*sinθ2 - cosθ2*sinθ1), cosθ4*(-cosθ1*sinθ2 - cosθ2*sinθ1) - sinθ4*(cosθ1*cosθ2 - sinθ1*sinθ2), 0,275*cosθ1*cosθ2 + 325*cosθ1 - 275*sinθ1*sinθ2]
 [cosθ4*(cosθ1*sinθ2 + cosθ2*sinθ1) + sinθ4*(cosθ1*cosθ2 - sinθ1*sinθ2), cosθ4*(cosθ1*cosθ2 - sinθ1*sinθ2) - sinθ4*(cosθ1*sinθ2 + cosθ2*sinθ1), 0, 275*cosθ1*sinθ2 + 275*cosθ2*sinθ1 + 325*sinθ1]
 [0, 0, 1, J3 + 150]
 [0,0,0,1]]
'''

'''
t = sympy.Matrix([[sinθ1θ2θ4,cosθ1θ2θ4-1,0],[cosθ1θ2θ4-1,sinθ1θ2θ4,0], [0,0,θ1θ2θ4]])
V = sympy.Matrix([[v1],[v2],[v3]])
P = sympy.Matrix([[px],[py],[pz]])


r2 = sympy.solve(t*V-P, [v1,v2,v3])
print('r2 = ', r2)
v1: (cosθ1θ2θ4*py - px*sinθ1θ2θ4 - py)/(cosθ1θ2θ4**2 - 2*cosθ1θ2θ4 - sinθ1θ2θ4**2 + 1)
v2: (cosθ1θ2θ4*px - px - py*sinθ1θ2θ4)/(cosθ1θ2θ4**2 - 2*cosθ1θ2θ4 - sinθ1θ2θ4**2 + 1)
v3: pz/θ1θ2θ4
'''


t = sympy.Matrix([[-ox,oy-1],[oy-1,-ox]])
V__2 = sympy.Matrix([[V1],[V2]])
p = sympy.Matrix([[Px],[Py]])
r2 = sympy.solve([-ox*V1+V2*(oy-1)-Px, (oy-1)*V1-V2*ox-Py ], [v1,v2])
print('r2 = ', r2)

'''MZ07转换矩阵
T06 =  [[cosθ6*(cosθ5*(cosθ4*(cosθ1*cosθ2*cosθ3 - cosθ1*sinθ2*sinθ3) + sinθ1*sinθ4) - sinθ5*(cosθ1*cosθ2*sinθ3 + cosθ1*cosθ3*sinθ2)) + sinθ6*(cosθ4*sinθ1 - sinθ4*(cosθ1*cosθ2*cosθ3 - cosθ1*sinθ2*sinθ3)), 
cosθ6*(cosθ4*sinθ1 - sinθ4*(cosθ1*cosθ2*cosθ3 - cosθ1*sinθ2*sinθ3)) - sinθ6*(cosθ5*(cosθ4*(cosθ1*cosθ2*cosθ3 - cosθ1*sinθ2*sinθ3) + sinθ1*sinθ4) - sinθ5*(cosθ1*cosθ2*sinθ3 + cosθ1*cosθ3*sinθ2)), 
cosθ5*(cosθ1*cosθ2*sinθ3 + cosθ1*cosθ3*sinθ2) + sinθ5*(cosθ4*(cosθ1*cosθ2*cosθ3 - cosθ1*sinθ2*sinθ3) + sinθ1*sinθ4), 
45*cosθ1*cosθ2*cosθ3 + 340*cosθ1*cosθ2*sinθ3 + 330*cosθ1*cosθ2 + 340*cosθ1*cosθ3*sinθ2 - 45*cosθ1*sinθ2*sinθ3 + 50*cosθ1 + 73*cosθ5*(cosθ1*cosθ2*sinθ3 + cosθ1*cosθ3*sinθ2) + 73*sinθ5*(cosθ4*(cosθ1*cosθ2*cosθ3 - cosθ1*sinθ2*sinθ3) + sinθ1*sinθ4)], 

[cosθ6*(cosθ5*(-cosθ1*sinθ4 + cosθ4*(cosθ2*cosθ3*sinθ1 - sinθ1*sinθ2*sinθ3)) - sinθ5*(cosθ2*sinθ1*sinθ3 + cosθ3*sinθ1*sinθ2)) + sinθ6*(-cosθ1*cosθ4 - sinθ4*(cosθ2*cosθ3*sinθ1 - sinθ1*sinθ2*sinθ3)), 
cosθ6*(-cosθ1*cosθ4 - sinθ4*(cosθ2*cosθ3*sinθ1 - sinθ1*sinθ2*sinθ3)) - sinθ6*(cosθ5*(-cosθ1*sinθ4 + cosθ4*(cosθ2*cosθ3*sinθ1 - sinθ1*sinθ2*sinθ3)) - sinθ5*(cosθ2*sinθ1*sinθ3 + cosθ3*sinθ1*sinθ2)), 
cosθ5*(cosθ2*sinθ1*sinθ3 + cosθ3*sinθ1*sinθ2) + sinθ5*(-cosθ1*sinθ4 + cosθ4*(cosθ2*cosθ3*sinθ1 - sinθ1*sinθ2*sinθ3)), 
45*cosθ2*cosθ3*sinθ1 + 340*cosθ2*sinθ1*sinθ3 + 330*cosθ2*sinθ1 + 340*cosθ3*sinθ1*sinθ2 + 73*cosθ5*(cosθ2*sinθ1*sinθ3 + cosθ3*sinθ1*sinθ2) - 45*sinθ1*sinθ2*sinθ3 + 50*sinθ1 + 73*sinθ5*(-cosθ1*sinθ4 + cosθ4*(cosθ2*cosθ3*sinθ1 - sinθ1*sinθ2*sinθ3))], 

[cosθ6*(cosθ4*cosθ5*(cosθ2*sinθ3 + cosθ3*sinθ2) - sinθ5*(-cosθ2*cosθ3 + sinθ2*sinθ3)) - sinθ4*sinθ6*(cosθ2*sinθ3 + cosθ3*sinθ2), 
-cosθ6*sinθ4*(cosθ2*sinθ3 + cosθ3*sinθ2) - sinθ6*(cosθ4*cosθ5*(cosθ2*sinθ3 + cosθ3*sinθ2) - sinθ5*(-cosθ2*cosθ3 + sinθ2*sinθ3)), 
cosθ4*sinθ5*(cosθ2*sinθ3 + cosθ3*sinθ2) + cosθ5*(-cosθ2*cosθ3 + sinθ2*sinθ3), 
-340*cosθ2*cosθ3 + 45*cosθ2*sinθ3 + 45*cosθ3*sinθ2 + 73*cosθ4*sinθ5*(cosθ2*sinθ3 + cosθ3*sinθ2) + 73*cosθ5*(-cosθ2*cosθ3 + sinθ2*sinθ3) + 340*sinθ2*sinθ3 + 330*sinθ2 + 345], 

[0, 0, 0, 1]]
'''
X = 45*cosθ1*cosθ2*cosθ3 + 340*cosθ1*cosθ2*sinθ3 + 330*cosθ1*cosθ2 + 340*cosθ1*cosθ3*sinθ2 - 45*cosθ1*sinθ2*sinθ3 + 50*cosθ1 + 73*cosθ5*(cosθ1*cosθ2*sinθ3 + cosθ1*cosθ3*sinθ2) + 73*sinθ5*(cosθ4*(cosθ1*cosθ2*cosθ3 - cosθ1*sinθ2*sinθ3) + sinθ1*sinθ4)

Y = 45*cosθ2*cosθ3*sinθ1 + 340*cosθ2*sinθ1*sinθ3 + 330*cosθ2*sinθ1 + 340*cosθ3*sinθ1*sinθ2 + 73*cosθ5*(cosθ2*sinθ1*sinθ3 + cosθ3*sinθ1*sinθ2) - 45*sinθ1*sinθ2*sinθ3 + 50*sinθ1 + 73*sinθ5*(-cosθ1*sinθ4 + cosθ4*(cosθ2*cosθ3*sinθ1 - sinθ1*sinθ2*sinθ3))

Z = -340*cosθ2*cosθ3 + 45*cosθ2*sinθ3 + 45*cosθ3*sinθ2 + 73*cosθ4*sinθ5*(cosθ2*sinθ3 + cosθ3*sinθ2) + 73*cosθ5*(-cosθ2*cosθ3 + sinθ2*sinθ3) + 340*sinθ2*sinθ3 + 330*sinθ2 + 345

e_