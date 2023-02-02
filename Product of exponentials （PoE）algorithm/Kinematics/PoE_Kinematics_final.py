
# -*- coding: utf-8 -*-
# coding=utf-8
# coding: utf-8

import sympy
import numpy as np
import matplotlib.pyplot as plt

"""RPY角旋转矩阵公式
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
"""
'''RPY角运算公式
RPY_R = Rz*Ry*Rx

RPY_Rotation =  sympy.Matrix([[cosθ*cosϕ, -cosψ*sinϕ + cosϕ*sinθ*sinψ, cosψ*cosϕ*sinθ + sinψ*sinϕ, 0], 
						 	 [cosθ*sinϕ, cosψ*cosϕ + sinθ*sinψ*sinϕ, cosψ*sinθ*sinϕ - cosϕ*sinψ, 0], 
						     [-sinθ, cosθ*sinψ, cosθ*cosψ, 0], 
						 	 [0, 0, 0, 1]])
'''
'''
ω  = sympy.Matrix([[0, -ω3, ω2], 
             	   [ω3, 0, -ω1],
			 	   [-ω2, ω1, 0]])

Rodrigues’s formula
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
		I = np.array([[1,0,0],
					 [0,1,0],
					 [0,0,1]])

		ω = np.array([[0,-ω3,ω2],
					 [ω3,0,-ω1],
					 [-ω2,ω1,0]])

		R = I + np.dot(np.sin(θ),ω) + np.dot((1-np.cos(θ)),np.dot(ω,ω))

		return(R)
	def e_Sθ(ω1,ω2,ω3,v1,v2,v3,θ):
		ω = np.array([[0,-ω3,ω2],
					 [ω3,0,-ω1],
					 [-ω2,ω1,0]])

		I = np.array([[1,0,0],
					 [0,1,0],
					 [0,0,1]])

		v = np.random.random(3)
		v = v.reshape(3, 1)
		v[0,0] = v1                     
		v[1,0] = v2                
		v[2,0] = v3              

		e_Sθ = np.random.random(16)
		e_Sθ = e_Sθ.reshape(4, 4)

		R = PoE.e_ωθ(ω1,ω2,ω3,θ)

		e_Sθ[0,0] = R[0,0]
		e_Sθ[1,0] = R[1,0]
		e_Sθ[2,0] = R[2,0]
		e_Sθ[3,0] = 0

		e_Sθ[0,1] = R[0,1]
		e_Sθ[1,1] = R[1,1]
		e_Sθ[2,1] = R[2,1]
		e_Sθ[3,1] = 0

		e_Sθ[0,2] = R[0,2]
		e_Sθ[1,2] = R[1,2]
		e_Sθ[2,2] = R[2,2]
		e_Sθ[3,2] = 0

		V = np.dot(np.dot(I,θ)+np.dot((1-np.cos(θ)),ω)+np.dot((θ-np.sin(θ)),np.dot(ω,ω)),v)

		e_Sθ[0,3] = V[0]
		e_Sθ[1,3] = V[1]
		e_Sθ[2,3] = V[2]
		e_Sθ[3,3] = 1
		
		return e_Sθ
	def AdT(ω1,ω2,ω3,v1,v2,v3,θ):
		T = PoE.e_Sθ(ω1,ω2,ω3,v1,v2,v3,θ)

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

	M = np.random.random(16)
	M = M.reshape(4, 4)
	M[0,0] = 1
	M[1,0] = 0
	M[2,0] = 0
	M[3,0] = 0

	M[0,1] = 0
	M[1,1] = -1
	M[2,1] = 0
	M[3,1] = 0

	M[0,2] = 0
	M[1,2] = 0
	M[2,2] = -1
	M[3,2] = 0

	M[0,3] = 50+330+45
	M[1,3] = 0
	M[2,3] = 345-340-73
	M[3,3] = 1

	S1 = np.array([[0],[0],[1],[0],[0],[0]])
	S2 = np.array([[0],[-1],[0],[345],[0],[-50]])
	S3 = np.array([[0],[-1],[0],[345],[0],[-330-50]])
	S4 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])
	S5 = np.array([[0],[-1],[0],[345-340],[0],[-50-330-45]])
	S6 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])
	
	e_Sθ1 = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ1*np.pi/180)
	e_Sθ2 = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ2*np.pi/180)
	e_Sθ3 = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ3*np.pi/180)
	e_Sθ4 = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ4*np.pi/180)
	e_Sθ5 = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],θ5*np.pi/180)
	e_Sθ6 = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],θ6*np.pi/180)

	T12 = np.dot(e_Sθ1, e_Sθ2)
	T13 = np.dot(T12, e_Sθ3)
	T14 = np.dot(T13, e_Sθ4)
	T15 = np.dot(T14, e_Sθ5)
	T16 = np.dot(T15, e_Sθ6)
	T = np.dot(T16, M) 

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

	'''
	求矩阵的逆
	Js_Inverse = np.linalg.inv(Js)
	求矩阵的伪逆(瘦高型)
	Js_T = Js.T	
	Js_pseudoInverse = np.dot(np.linalg.inv(np.dot(Js_T,Js)),Js_T)
	'''
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
def MZ07_Forward_Kinematics_PoEModel_EndEffectorFrame(θ1,θ2,θ3,θ4,θ5,θ6):
	M = np.array([[1,0,0,50+330+45],
				  [0,-1,0,0],
				  [0,0,-1,345-340-73],
				  [0,0,0,1]])

	B1 = np.array([[0],[0],[-1],[0],[-50-330-45],[0]])
	B2 = np.array([[0],[1],[0],[340+73],[0],[-330-45]])
	B3 = np.array([[0],[1],[0],[340+73],[0],[-45]])
	B4 = np.array([[0],[0],[1],[0],[0],[0]])
	B5 = np.array([[0],[1],[0],[73],[0],[0]])	
	B6 = np.array([[0],[0],[1],[0],[0],[0]])

	e_Bθ1 = PoE.e_Sθ(B1[0][0],B1[1][0],B1[2][0],B1[3][0],B1[4][0],B1[5][0],θ1*np.pi/180)
	e_Bθ2 = PoE.e_Sθ(B2[0][0],B2[1][0],B2[2][0],B2[3][0],B2[4][0],B2[5][0],θ2*np.pi/180)
	e_Bθ3 = PoE.e_Sθ(B3[0][0],B3[1][0],B3[2][0],B3[3][0],B3[4][0],B3[5][0],θ3*np.pi/180)
	e_Bθ4 = PoE.e_Sθ(B4[0][0],B4[1][0],B4[2][0],B4[3][0],B4[4][0],B4[5][0],θ4*np.pi/180)
	e_Bθ5 = PoE.e_Sθ(B5[0][0],B5[1][0],B5[2][0],B5[3][0],B5[4][0],B5[5][0],θ5*np.pi/180)
	e_Bθ6 = PoE.e_Sθ(B6[0][0],B6[1][0],B6[2][0],B6[3][0],B6[4][0],B6[5][0],θ6*np.pi/180)

	T12 = np.dot(e_Bθ1, e_Bθ2)
	T13 = np.dot(T12, e_Bθ3)
	T14 = np.dot(T13, e_Bθ4)
	T15 = np.dot(T14, e_Bθ5)
	T16 = np.dot(T15, e_Bθ6)
	T = np.dot(M,T16)

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
	print('T', T)
	return [px, py, pz, RPY_phi*180/np.pi, RPY_theta*180/np.pi, RPY_psi*180/np.pi, Euler_phi*180/np.pi,Euler_theta*180/np.pi,Euler_psi*180/np.pi]
def MZ07_Forward_Velocity_Kinematics_PoEModel(theta1,theta2,theta3,theta4,theta5,theta6,omega_θ1,omega_θ2,omega_θ3,omega_θ4,omega_θ5,omega_θ6):

	M = np.random.random(16)
	M = M.reshape(4, 4)
	M[0,0] = 1
	M[1,0] = 0
	M[2,0] = 0
	M[3,0] = 0

	M[0,1] = 0
	M[1,1] = -1
	M[2,1] = 0
	M[3,1] = 0

	M[0,2] = 0
	M[1,2] = 0
	M[2,2] = -1
	M[3,2] = 0

	M[0,3] = 50+330+45
	M[1,3] = 0
	M[2,3] = 345-340-73
	M[3,3] = 1

	S1 = np.array([[0],[0],[1],[0],[0],[0]])
	S2 = np.array([[0],[-1],[0],[345],[0],[-50]])
	S3 = np.array([[0],[-1],[0],[345],[0],[-330-50]])
	S4 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])
	S5 = np.array([[0],[-1],[0],[345-340],[0],[-50-330-45]])
	S6 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])

	AdT1 = PoE.AdT(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],theta1*np.pi/180)
	AdT2 = PoE.AdT(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],theta2*np.pi/180)
	AdT3 = PoE.AdT(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],theta3*np.pi/180)
	AdT4 = PoE.AdT(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],theta4*np.pi/180)
	AdT5 = PoE.AdT(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],theta5*np.pi/180)
	AdT6 = PoE.AdT(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],theta6*np.pi/180)

	AdT12 = np.dot(AdT1,AdT2)
	AdT123 = np.dot(AdT12,AdT3)
	AdT1234 = np.dot(AdT123,AdT4)
	AdT12345 = np.dot(AdT1234,AdT5)

	Js1 = S1
	Js2 = np.dot(AdT1,S2)
	Js3 = np.dot(AdT12,S3)
	Js4 = np.dot(AdT123,S4)
	Js5 = np.dot(AdT1234,S5)
	Js6 = np.dot(AdT12345,S6)

	Js = np.random.random(36)
	Js = Js.reshape(6, 6)

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

	Js[0,2] = Js3[0][0]
	Js[1,2] = Js3[1][0]
	Js[2,2] = Js3[2][0]
	Js[3,2] = Js3[3][0]
	Js[4,2] = Js3[4][0]
	Js[5,2] = Js3[5][0]

	Js[0,3] = Js4[0][0]
	Js[1,3] = Js4[1][0]
	Js[2,3] = Js4[2][0]
	Js[3,3] = Js4[3][0]
	Js[4,3] = Js4[4][0]
	Js[5,3] = Js4[5][0]

	Js[0,4] = Js5[0][0]
	Js[1,4] = Js5[1][0]
	Js[2,4] = Js5[2][0]
	Js[3,4] = Js5[3][0]
	Js[4,4] = Js5[4][0]
	Js[5,4] = Js5[5][0]

	Js[0,5] = Js6[0][0]
	Js[1,5] = Js6[1][0]
	Js[2,5] = Js6[2][0]
	Js[3,5] = Js6[3][0]
	Js[4,5] = Js6[4][0]
	Js[5,5] = Js6[5][0]

	omega_θ = [[omega_θ1],[omega_θ2],[omega_θ3],[omega_θ4],[omega_θ5],[omega_θ6]]

	V = np.dot(Js,omega_θ)
	print('角速度 =', omega_θ)
	print('运动旋量 = ', V)

	return [V[0][0],V[1][0],V[2][0],V[3][0],V[4][0],V[5][0]]

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

	'''
	e_Sθ1 = PoE.e_Sθ(0,0,1,0,0,0,θ1*np.pi/180)
	e_Sθ2 = PoE.e_Sθ(0,0,1,0,-325,0,θ2*np.pi/180)
	e_Sθ3 = PoE.e_Sθ(0,0,0,0,0,1,θ3)
	e_Sθ4 = PoE.e_Sθ(0,0,1,0,-325-275,0,θ4*np.pi/180)
	'''

	S1 = np.array([[0],[0],[1],[0],[0],[0]])
	S2 = np.array([[0],[0],[1],[0],[-325],[0]])
	S3 = np.array([[0],[0],[0],[0],[0],[1]])
	S4 = np.array([[0],[0],[1],[0],[-325-275],[0]])

	S1_T = S1.T[0]
	S2_T = S2.T[0]
	S3_T = S3.T[0]
	S4_T = S4.T[0]
	
	e_Sθ1 = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ1*np.pi/180)
	e_Sθ2 = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ2*np.pi/180)
	e_Sθ3 = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ3)
	e_Sθ4 = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ4*np.pi/180)

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
	
	Euler_phi = RPY_phi = np.arctan(ny/nx)

	return [px, py, pz, RPY_phi*180/np.pi, Euler_phi*180/np.pi]
def EC06_Forward_Kinematics_PoEModel_EndEffectorFrame(θ1,θ2,θ3,θ4):
	M = np.array([[1,0,0,325+275],
				  [0,1,0,0],
				  [0,0,1,150],
				  [0,0,0,1]])
	
	B1 = np.array([[0],[0],[1],[0],[325+275],[0]])
	B2 = np.array([[0],[0],[1],[0],[275],[0]])
	B3 = np.array([[0],[0],[0],[0],[0],[1]])
	B4 = np.array([[0],[0],[1],[0],[0],[0]])
	
	e_Bθ1 = PoE.e_Sθ(B1[0][0],B1[1][0],B1[2][0],B1[3][0],B1[4][0],B1[5][0],θ1*np.pi/180)
	e_Bθ2 = PoE.e_Sθ(B2[0][0],B2[1][0],B2[2][0],B2[3][0],B2[4][0],B2[5][0],θ2*np.pi/180)
	e_Bθ3 = PoE.e_Sθ(B3[0][0],B3[1][0],B3[2][0],B3[3][0],B3[4][0],B3[5][0],θ3)
	e_Bθ4 = PoE.e_Sθ(B4[0][0],B4[1][0],B4[2][0],B4[3][0],B4[4][0],B4[5][0],θ4*np.pi/180)

	T12 = np.dot(e_Bθ1, e_Bθ2)
	T13 = np.dot(T12, e_Bθ3)
	T14 = np.dot(T13, e_Bθ4)	
	T = np.dot(M,T14) 

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
	
	Euler_phi = RPY_phi = np.arctan(ny/nx)

	return [px, py, pz, RPY_phi*180/np.pi, Euler_phi*180/np.pi]
def EC06_Forward_Velocity_Kinematics_PoEModel(theta1,theta2,theta3,theta4,omega_θ1,omega_θ2,omega_θ3,omega_θ4):

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

	S1 = np.array([[0],[0],[1],[0],[0],[0]])
	S2 = np.array([[0],[0],[1],[0],[-325],[0]])
	S3 = np.array([[0],[0],[0],[0],[0],[1]])
	S4 = np.array([[0],[0],[1],[0],[-325-275],[0]])

	theta = np.array([[theta1*np.pi/180],
					[theta2*np.pi/180],
					[theta3],
					[theta4*np.pi/180]])

	AdT1 = PoE.AdT(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],theta[0][0])
	AdT2 = PoE.AdT(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],theta[1][0])
	AdT3 = PoE.AdT(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],theta[2][0])

	AdT12 = np.dot(AdT1,AdT2)
	AdT123 = np.dot(AdT12,AdT3)

	Js1 = S1
	Js2 = np.dot(AdT1,S2)
	Js3 = np.dot(AdT12,S3)
	Js4 = np.dot(AdT123,S4)

	Js = np.random.random(24)
	Js = Js.reshape(6, 4)

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

	Js[0,2] = Js3[0][0]
	Js[1,2] = Js3[1][0]
	Js[2,2] = Js3[2][0]
	Js[3,2] = Js3[3][0]
	Js[4,2] = Js3[4][0]
	Js[5,2] = Js3[5][0]

	Js[0,3] = Js4[0][0]
	Js[1,3] = Js4[1][0]
	Js[2,3] = Js4[2][0]
	Js[3,3] = Js4[3][0]
	Js[4,3] = Js4[4][0]
	Js[5,3] = Js4[5][0]

	omega_θ = [[omega_θ1],
			[omega_θ2],
			[omega_θ3],
			[omega_θ4]]

	'''
	J = np.array([[0,0,0,0],
				[0,0,0,0],
				[1,1,0,1],
				[0,325*np.sin(theta[0][0]),0, 325*np.sin(theta[0][0])+275*np.sin(theta[0][0]+theta[1][0])],
				[0,-325*np.cos(theta[0][0]),0,-325*np.cos(theta[0][0])-275*np.cos(theta[0][0]+theta[1][0])],
				[0,0,1,0]])
	'''

	V = np.dot(Js,omega_θ)
	print('角速度 =', omega_θ)
	print('运动旋量 = ', V)

	return [V[0][0],V[1][0],V[2][0],V[3][0]]
def EC06_Numerical_Inverse_Kinematics_PoEModel(X,Y,Z,θ,steps=10,error=0.001, theta0_1=5,theta0_2=5,theta0_3=150,theta0_4=5):

	θ = θ*np.pi/180

	'''目标位置'''
	Ts_target = np.array([[np.cos(θ), -np.sin(θ),  0, X], 
            	  		  [np.sin(θ),  np.cos(θ),  0 ,Y],
				 		  [0,           0 ,        1, Z],
				 		  [0,           0,         0 ,1]])

	M = np.array([[1,0,0,325+275],
				  [0,1,0,      0],
				  [0,0,1,    150],
				  [0,0,0,      1]])

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

	S1 = np.array([[0],[0],[1],[0],[0],[0]])
	S2 = np.array([[0],[0],[1],[0],[-325],[0]])
	S3 = np.array([[0],[0],[0],[0],[0],[1]])
	S4 = np.array([[0],[0],[1],[0],[-325-275],[0]])
	
	B1 = np.array([[0],[0],[1],[0],[325+275],[0]])
	B2 = np.array([[0],[0],[1],[0],[275],[0]])
	B3 = np.array([[0],[0],[0],[0],[0],[1]])
	B4 = np.array([[0],[0],[1],[0],[0],[0]])

	errorlist = []
	theta1 = []
	theta2 = []
	theta3 = []
	theta4 = []

	theta = np.array([[theta0_1*np.pi/180],
					  [theta0_2*np.pi/180],
					  [theta0_3],
					  [theta0_4*np.pi/180]])

	'''
	T04 =  [[cosθ4*(cosθ1*cosθ2 - sinθ1*sinθ2) + sinθ4*(-cosθ1*sinθ2 - cosθ2*sinθ1), cosθ4*(-cosθ1*sinθ2 - cosθ2*sinθ1) - sinθ4*(cosθ1*cosθ2 - sinθ1*sinθ2), 0,275*cosθ1*cosθ2 + 325*cosθ1 - 275*sinθ1*sinθ2]
	[cosθ4*(cosθ1*sinθ2 + cosθ2*sinθ1) + sinθ4*(cosθ1*cosθ2 - sinθ1*sinθ2), cosθ4*(cosθ1*cosθ2 - sinθ1*sinθ2) - sinθ4*(cosθ1*sinθ2 + cosθ2*sinθ1), 0, 275*cosθ1*sinθ2 + 275*cosθ2*sinθ1 + 325*sinθ1]
	[0, 0, 1, J3 + 150]
	[0,0,0,1]]
	'''

	'''
	t = sympy.Matrix([[sinθ1θ2θ4,cosθ1θ2θ4-1,0],[cosθ1θ2θ4-1,sinθ1θ2θ4,0], [0,0,θ1θ2θ4]])
	V = sympy.Matrix([[v1],[v2],[v3]])
	P = sympy.Matrix([[X],[Y],[Z]])

	r2 = sympy.solve(t*V-P, [v1,v2,v3])
	print('r2 = ', r2)
	v1: (cosθ1θ2θ4*Y - X*sinθ1θ2θ4 - Y)/(cosθ1θ2θ4**2 - 2*cosθ1θ2θ4 - sinθ1θ2θ4**2 + 1)
	v2: (cosθ1θ2θ4*X - X - Y*sinθ1θ2θ4)/(cosθ1θ2θ4**2 - 2*cosθ1θ2θ4 - sinθ1θ2θ4**2 + 1)
	v3: Z/θ1θ2θ4
	'''
	fw1 = open(r"C:\Users\eastw\Desktop\development\Product of exponentials （PoE）algorithm\1.txt", 'w')
	fw2 = open(r"C:\Users\eastw\Desktop\development\Product of exponentials （PoE）algorithm\2.txt", 'w')
	fw3 = open(r"C:\Users\eastw\Desktop\development\Product of exponentials （PoE）algorithm\3.txt", 'w')
	fw4 = open(r"C:\Users\eastw\Desktop\development\Product of exponentials （PoE）algorithm\4.txt", 'w')

	for i in range(steps): 

		e_Sθ1 = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],theta[0][0])
		e_Sθ2 = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],theta[1][0])
		e_Sθ3 = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],theta[2][0])
		e_Sθ4 = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],theta[3][0])

		Ts_b12 = np.dot(e_Sθ1, e_Sθ2)
		Ts_b13 = np.dot(Ts_b12, e_Sθ3)
		Ts_b14 = np.dot(Ts_b13, e_Sθ4)	
		Ts_b = np.dot(Ts_b14,M) 

		Ts_b_Inverse = np.linalg.inv(Ts_b)		
		Tb_s = Ts_b_Inverse 
		Tb_target = np.dot(Tb_s,Ts_target)

		nx = Tb_target[0,0]
		ny = Tb_target[1,0]
		nz = Tb_target[2,0]
		ox = Tb_target[0,1]
		oy = Tb_target[1,1]
		oz = Tb_target[2,1]
		ax = Tb_target[0,2]
		ay = Tb_target[1,2]
		az = Tb_target[2,2]
		px = Tb_target[0,3]
		py = Tb_target[1,3]
		pz = Tb_target[2,3]

		θi = np. arccos((nx + oy + az  - 1)/2)

		ωi1 = (oz-ay)/(2*np.sin(θi))
		ωi2 = (ax-nz)/(2*np.sin(θi))
		ωi3 = (ny-ox)/(2*np.sin(θi))

		print('ωi3 = ',ωi3 )
		ωi = np.random.random(9)
		ωi = ωi.reshape(3, 3)

		ωi[0,0] = 0               
		ωi[0,1] = -ωi3                   
		ωi[0,2] = ωi2        
		ωi[1,0] = ωi3              
		ωi[1,1] = 0                  
		ωi[1,2] = -ωi1    
		ωi[2,0] = -ωi2               
		ωi[2,1] = ωi1                   
		ωi[2,2] = 0  


		Ri = np.dot(I,θi)+np.dot((1-np.cos(θi)),ωi)+np.dot((θi-np.sin(θi)),np.dot(ωi,ωi))		
		print('Ri = ',Ri)	
		Ri_Inverse = np.linalg.inv(Ri)

		Pi = [[px],
			  [py],
			  [pz]]
		vi = np.dot(Ri_Inverse,Pi)

		vi1 = vi[0][0]
		vi2 = vi[1][0]
		vi3 = vi[2][0]

		Vbi = np.array([[0],
						[0],
						[θi],
						[vi1*θi],
						[vi2*θi],
						[vi3*θi]])

		AdT_b1 = PoE.AdT(B1[0][0],B1[1][0],B1[2][0],B1[3][0],B1[4][0],B1[5][0],-theta[0][0])
		AdT_b2 = PoE.AdT(B2[0][0],B2[1][0],B2[2][0],B2[3][0],B2[4][0],B2[5][0],-theta[1][0])
		AdT_b3 = PoE.AdT(B3[0][0],B3[1][0],B3[2][0],B3[3][0],B3[4][0],B3[5][0],-theta[2][0])
		AdT_b4 = PoE.AdT(B4[0][0],B4[1][0],B4[2][0],B4[3][0],B4[4][0],B4[5][0],-theta[3][0])

		AdT_b43 = np.dot(AdT_b4,AdT_b3)
		AdT_b432 = np.dot(AdT_b43,AdT_b2)

		Jb1 = np.dot(AdT_b432,B1)
		Jb2 = np.dot(AdT_b43,B2)
		Jb3 = np.dot(AdT_b4,B3)
		Jb4 = B4

		Jb = np.random.random(24)
		Jb = Jb.reshape(6, 4)

		Jb[0,0] = Jb1[0][0]
		Jb[1,0] = Jb1[1][0]
		Jb[2,0] = Jb1[2][0]
		Jb[3,0] = Jb1[3][0]
		Jb[4,0] = Jb1[4][0]
		Jb[5,0] = Jb1[5][0]

		Jb[0,1] = Jb2[0][0]
		Jb[1,1] = Jb2[1][0]
		Jb[2,1] = Jb2[2][0]
		Jb[3,1] = Jb2[3][0]
		Jb[4,1] = Jb2[4][0]
		Jb[5,1] = Jb2[5][0]

		Jb[0,2] = Jb3[0][0]
		Jb[1,2] = Jb3[1][0]
		Jb[2,2] = Jb3[2][0]
		Jb[3,2] = Jb3[3][0]
		Jb[4,2] = Jb3[4][0]
		Jb[5,2] = Jb3[5][0]

		Jb[0,3] = Jb4[0][0]
		Jb[1,3] = Jb4[1][0]
		Jb[2,3] = Jb4[2][0]
		Jb[3,3] = Jb4[3][0]
		Jb[4,3] = Jb4[4][0]
		Jb[5,3] = Jb4[5][0]

		'''求矩阵的伪逆(瘦高型)'''	
		Jb_T = Jb.T	
		Js_pseudoInverse = np.dot(np.linalg.inv(np.dot(Jb_T,Jb)),Jb_T)
		thetaX = theta
		theta = theta + np.dot(Js_pseudoInverse,Vbi)
		deltatheta = theta-thetaX

		e = np.sqrt(np.dot(deltatheta[0][0],deltatheta[0][0]) + np.dot(deltatheta[1][0] ,deltatheta[1][0] ) + np.dot(deltatheta[2][0] ,deltatheta[2][0] ) + np.dot(deltatheta[3][0] ,deltatheta[3][0]))
	
		'''
		e124 = np.sqrt(np.dot(deltatheta[0][0],deltatheta[0][0]) + np.dot(deltatheta[1][0],deltatheta[1][0]) + np.dot(deltatheta[3][0],deltatheta[3][0]))
		e3 = deltatheta[3][0]
		e = e124 + e3
		'''
		errorlist.append(e)

		theta1.append(theta[0][0]*180/np.pi)
		theta2.append(theta[1][0]*180/np.pi)
		theta3.append(theta[2][0])
		theta4.append(theta[3][0]*180/np.pi)

		str1 = str(theta[0][0]*180/np.pi)

		fw1.write(str1)
		fw1.write('\n')
		str2 = str(theta[1][0]*180/np.pi)

		fw2.write(str2)
		fw2.write('\n')
		str3 = str(theta[2][0])
	
		fw3.write(str3)
		fw3.write('\n')
		str4 = str(theta[3][0]*180/np.pi)
		
		fw4.write(str4)
		fw4.write('\n')
		
		if errorlist[-1] < error:
			print('iteration: ', theta)
			
			break
	x_axis = [i for i in range(len(errorlist))]
	
	plt.figure()
	plta = plt.gca(projection='3d')
	plta.plot(x_axis, errorlist,0)
	plta.plot(x_axis, theta1,1)	
	plta.plot(x_axis, theta2,2)
	plta.plot(x_axis, theta3,3)
	plta.plot(x_axis, theta4,4)

	plt.show()
	
	return [theta[0][0]*180/np.pi,theta[1][0]*180/np.pi,theta[2][0],theta[3][0]*180/np.pi]
def EC06_Trajectory_Planning_StraightLine_Polynomial_AngleOutput(θ1_origin, θ2_origin,θ3_origin,θ4_origin,θ1_target, θ2_target,θ3_target,θ4_target,T=3,t=0.01,steps=10,error=0.001):


	M = np.array([[1,0,0,325+275],
				  [0,1,0,      0],
				  [0,0,1,    150],
				  [0,0,0,      1]])
	
	errorlist = []
	theta1 = []
	theta2 = []
	theta3 = []
	theta4 = []

	I = np.array([[1,0,0],
				  [0,1,0],
				  [0,0,1]])

	S1 = np.array([[0],[0],[1],[0],[0],[0]])
	S2 = np.array([[0],[0],[1],[0],[-325],[0]])
	S3 = np.array([[0],[0],[0],[0],[0],[1]])
	S4 = np.array([[0],[0],[1],[0],[-325-275],[0]])
	
	B1 = np.array([[0],[0],[1],[0],[325+275],[0]])
	B2 = np.array([[0],[0],[1],[0],[275],[0]])
	B3 = np.array([[0],[0],[0],[0],[0],[1]])
	B4 = np.array([[0],[0],[1],[0],[0],[0]])

	errorlist = []
	theta1 = []
	theta2 = []
	theta3 = []
	theta4 = []

	'''设定终点'''
	θ_target = np.array([[θ1_target*np.pi/180],
					  	 [θ2_target*np.pi/180],
						 [θ3_target],
						 [θ4_target*np.pi/180]])

	e_Sθ1_target = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ_target[0][0])
	e_Sθ2_target = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ_target[1][0])
	e_Sθ3_target = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ_target[2][0])
	e_Sθ4_target = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ_target[3][0])


	T12_target = np.dot(e_Sθ1_target, e_Sθ2_target)
	T13_target = np.dot(T12_target, e_Sθ3_target)
	T14_target = np.dot(T13_target, e_Sθ4_target)
	Ts_target = np.dot(T14_target, M) 

	'''设定起点'''
	θ_origin = np.array([[θ1_origin*np.pi/180],
						 [θ2_origin*np.pi/180],
						 [θ3_origin],
						 [θ4_origin*np.pi/180]])

	e_Sθ1_origin = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ_origin[0][0])
	e_Sθ2_origin = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ_origin[1][0])
	e_Sθ3_origin = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ_origin[2][0])
	e_Sθ4_origin = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ_origin[3][0])

	T12_origin = np.dot(e_Sθ1_origin, e_Sθ2_origin)
	T13_origin = np.dot(T12_origin, e_Sθ3_origin)
	T14_origin = np.dot(T13_origin, e_Sθ4_origin)
	Ts_origin = np.dot(T14_origin, M) 

	Px_Path = []
	Py_Path = []
	Pz_Path = []
	vx_Path = []
	vy_Path = []
	vz_Path = []
	v_Path = []
	for i in np.arange(0,T+t,t):
		if i == 0:

			theta = θ_origin
			theta1.append(theta[0][0]*180/np.pi)
			theta2.append(theta[1][0]*180/np.pi)
			theta3.append(theta[2][0])
			theta4.append(theta[3][0]*180/np.pi)


			Px_Path.append(Ts_origin[0,3])
			Py_Path.append(Ts_origin[1,3])
			Pz_Path.append(Ts_origin[2,3])

			V_Tst = (6*i/(T*T)-6*i*i/(T*T*T))*(Ts_target-Ts_origin)

			vx_Path.append(V_Tst[0,3])
			vy_Path.append(V_Tst[1,3])
			vz_Path.append(V_Tst[2,3])

			v_Path.append(np.sqrt(vx_Path[-1]*vx_Path[-1] + vy_Path[-1]*vy_Path[-1] + vz_Path[-1]*vz_Path[-1]))
			
			i = i + t

		Tst = Ts_origin + (3*i*i/(T*T)-2*i*i*i/(T*T*T))*(Ts_target-Ts_origin)

		'''逆运动学数值解'''
		for j in range(steps): 

			e_Sθ1 = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],theta[0][0])
			e_Sθ2 = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],theta[1][0])
			e_Sθ3 = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],theta[2][0])
			e_Sθ4 = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],theta[3][0])

			Ts_b12 = np.dot(e_Sθ1, e_Sθ2)
			Ts_b13 = np.dot(Ts_b12, e_Sθ3)
			Ts_b14 = np.dot(Ts_b13, e_Sθ4)	
			Ts_b = np.dot(Ts_b14,M) 

			Ts_b_Inverse = np.linalg.inv(Ts_b)		
			Tb_s = Ts_b_Inverse 
			Tb_target = np.dot(Tb_s,Tst)

			nx = Tb_target[0,0]
			ny = Tb_target[1,0]
			nz = Tb_target[2,0]
			ox = Tb_target[0,1]
			oy = Tb_target[1,1]
			oz = Tb_target[2,1]
			ax = Tb_target[0,2]
			ay = Tb_target[1,2]
			az = Tb_target[2,2]
			px = Tb_target[0,3]
			py = Tb_target[1,3]
			pz = Tb_target[2,3]

			θi = np. arccos((nx + oy + az  - 1)/2)


			ωi1 = (oz-ay)/(2*np.sin(θi))
			ωi2 = (ax-nz)/(2*np.sin(θi))
			ωi3 = (ny-ox)/(2*np.sin(θi))

			ωi = np.random.random(9)
			ωi = ωi.reshape(3, 3)

			ωi[0,0] = 0               
			ωi[0,1] = -ωi3                   
			ωi[0,2] = ωi2        
			ωi[1,0] = ωi3              
			ωi[1,1] = 0                  
			ωi[1,2] = -ωi1    
			ωi[2,0] = -ωi2               
			ωi[2,1] = ωi1                   
			ωi[2,2] = 0  

			'''
			Ri_Inverse = I/θi-ωi/2+ np.dot((1/θi-1/(np.tan(θi/2)*2)),np.dot(ωi,ωi))

			'''
			Ri = np.dot(I,θi)+np.dot((1-np.cos(θi)),ωi)+np.dot((θi-np.sin(θi)),np.dot(ωi,ωi))			

			Ri_Inverse = np.linalg.inv(Ri)

			Pi = [[px],
				  [py],
				  [pz]]
			vi = np.dot(Ri_Inverse,Pi)

			vi1 = vi[0][0]
			vi2 = vi[1][0]
			vi3 = vi[2][0]

			Vbi = np.array([[ωi1*θi],
							[ωi2*θi],
							[ωi3*θi],
							[vi1*θi],
							[vi2*θi],
							[vi3*θi]])
		
			AdT_b1 = PoE.AdT(B1[0][0],B1[1][0],B1[2][0],B1[3][0],B1[4][0],B1[5][0],-theta[0][0])
			AdT_b2 = PoE.AdT(B2[0][0],B2[1][0],B2[2][0],B2[3][0],B2[4][0],B2[5][0],-theta[1][0])
			AdT_b3 = PoE.AdT(B3[0][0],B3[1][0],B3[2][0],B3[3][0],B3[4][0],B3[5][0],-theta[2][0])
			AdT_b4 = PoE.AdT(B4[0][0],B4[1][0],B4[2][0],B4[3][0],B4[4][0],B4[5][0],-theta[3][0])

			AdT_b43 = np.dot(AdT_b4,AdT_b3)
			AdT_b432 = np.dot(AdT_b43,AdT_b2)

			Jb1 = np.dot(AdT_b432,B1)
			Jb2 = np.dot(AdT_b43,B2)
			Jb3 = np.dot(AdT_b4,B3)
			Jb4 = B4

			Jb = np.random.random(24)
			Jb = Jb.reshape(6, 4)

			Jb[0,0] = Jb1[0][0]
			Jb[1,0] = Jb1[1][0]
			Jb[2,0] = Jb1[2][0]
			Jb[3,0] = Jb1[3][0]
			Jb[4,0] = Jb1[4][0]
			Jb[5,0] = Jb1[5][0]

			Jb[0,1] = Jb2[0][0]
			Jb[1,1] = Jb2[1][0]
			Jb[2,1] = Jb2[2][0]
			Jb[3,1] = Jb2[3][0]
			Jb[4,1] = Jb2[4][0]
			Jb[5,1] = Jb2[5][0]

			Jb[0,2] = Jb3[0][0]
			Jb[1,2] = Jb3[1][0]
			Jb[2,2] = Jb3[2][0]
			Jb[3,2] = Jb3[3][0]
			Jb[4,2] = Jb3[4][0]
			Jb[5,2] = Jb3[5][0]

			Jb[0,3] = Jb4[0][0]
			Jb[1,3] = Jb4[1][0]
			Jb[2,3] = Jb4[2][0]
			Jb[3,3] = Jb4[3][0]
			Jb[4,3] = Jb4[4][0]
			Jb[5,3] = Jb4[5][0]

			'''求矩阵的伪逆(瘦高型)'''	
			Jb_T = Jb.T	
			Jb_pseudoInverse = np.dot(np.linalg.inv(np.dot(Jb_T,Jb)),Jb_T)
			thetaX = theta
			theta = theta + np.dot(Jb_pseudoInverse,Vbi)
			deltatheta = theta-thetaX

			#errorlist.append(np.linalg.norm(deltatheta))
			e = np.sqrt(np.dot(deltatheta[0][0],deltatheta[0][0]) + np.dot(deltatheta[1][0],deltatheta[1][0]) + np.dot(deltatheta[2][0],deltatheta[2][0]) + np.dot(deltatheta[3][0],deltatheta[3][0]))
			errorlist.append(e)
				
			if errorlist[-1] < error:
				# print('iteration: ', theta*180/np.pi%360)
				print('迭代次数:',j+1)
				theta1.append(theta[0][0]*180/np.pi)
				theta2.append(theta[1][0]*180/np.pi)
				theta3.append(theta[2][0])
				theta4.append(theta[3][0]*180/np.pi)
				break


		Px_Path.append(Tst[0,3])
		Py_Path.append(Tst[1,3])
		Pz_Path.append(Tst[2,3])


		V_Tst = (6*i/(T*T)-6*i*i/(T*T*T))*(Ts_target-Ts_origin)

		vx_Path.append(V_Tst[0,3])
		vy_Path.append(V_Tst[1,3])
		vz_Path.append(V_Tst[2,3])

		v_Path.append(np.sqrt(vx_Path[-1]*vx_Path[-1] + vy_Path[-1]*vy_Path[-1] + vz_Path[-1]*vz_Path[-1]))


	
	x_axis = [i for i in range(len(theta1))]
	x1_axis = [i for i in range(len(v_Path))]	
	plt.figure()
	# plta = plt.gca(projection='3d')
	# plta.plot(Px_Path, Py_Path, Pz_Path)
	plt.plot(x_axis,theta1)
	plt.plot(x_axis,theta2)
	plt.plot(x_axis,theta3)
	plt.plot(x_axis,theta4)

	# plt.plot(x1_axis,v_Path)

	plt.show()

	return np.array([[theta1],
					[theta2],
					[theta3],
					[theta4],
					[Px_Path],
					[Py_Path],
					[Pz_Path],
					[v_Path]])
def MZ07_Numerical_Inverse_Kinematics_PoEModel(X,Y,Z,Roll,Pitch,Yaw, steps = 10,error = 0.01,theta0_1=10,theta0_2=10,theta0_3=10,theta0_4=10,theta0_5=10,theta0_6=10):

	'''T06 =  Matrix([[cosθ6*(cosθ5*(cosθ4*(cosθ1*cosθ2*cosθ3 - cosθ1*sinθ2*sinθ3) + sinθ1*sinθ4) - sinθ5*(cosθ1*cosθ2*sinθ3 + cosθ1*cosθ3*sinθ2)) + sinθ6*(cosθ4*sinθ1 - sinθ4*(cosθ1*cosθ2*cosθ3 - cosθ1*sinθ2*sinθ3)), cosθ6*(cosθ4*sinθ1 - sinθ4*(cosθ1*cosθ2*cosθ3 - cosθ1*sinθ2*sinθ3)) - sinθ6*(cosθ5*(cosθ4*(cosθ1*cosθ2*cosθ3 - cosθ1*sinθ2*sinθ3) + sinθ1*sinθ4) - sinθ5*(cosθ1*cosθ2*sinθ3 + cosθ1*cosθ3*sinθ2)), cosθ5*(cosθ1*cosθ2*sinθ3 + cosθ1*cosθ3*sinθ2) + sinθ5*(cosθ4*(cosθ1*cosθ2*cosθ3 - cosθ1*sinθ2*sinθ3) + sinθ1*sinθ4), 45*cosθ1*cosθ2*cosθ3 + 340*cosθ1*cosθ2*sinθ3 + 330*cosθ1*cosθ2 + 340*cosθ1*cosθ3*sinθ2 - 45*cosθ1*sinθ2*sinθ3 + 50*cosθ1 + 73*cosθ5*(cosθ1*cosθ2*sinθ3 + cosθ1*cosθ3*sinθ2) + 73*sinθ5*(cosθ4*(cosθ1*cosθ2*cosθ3 - cosθ1*sinθ2*sinθ3) + sinθ1*sinθ4)], 
	[cosθ6*(cosθ5*(-cosθ1*sinθ4 + cosθ4*(cosθ2*cosθ3*sinθ1 - sinθ1*sinθ2*sinθ3)) - sinθ5*(cosθ2*sinθ1*sinθ3 + cosθ3*sinθ1*sinθ2)) + sinθ6*(-cosθ1*cosθ4 - sinθ4*(cosθ2*cosθ3*sinθ1 - sinθ1*sinθ2*sinθ3)), cosθ6*(-cosθ1*cosθ4 - sinθ4*(cosθ2*cosθ3*sinθ1 - sinθ1*sinθ2*sinθ3)) - sinθ6*(cosθ5*(-cosθ1*sinθ4 + cosθ4*(cosθ2*cosθ3*sinθ1 - sinθ1*sinθ2*sinθ3)) - sinθ5*(cosθ2*sinθ1*sinθ3 + cosθ3*sinθ1*sinθ2)), cosθ5*(cosθ2*sinθ1*sinθ3 + cosθ3*sinθ1*sinθ2) + sinθ5*(-cosθ1*sinθ4 + cosθ4*(cosθ2*cosθ3*sinθ1 - sinθ1*sinθ2*sinθ3)), 45*cosθ2*cosθ3*sinθ1 + 340*cosθ2*sinθ1*sinθ3 + 330*cosθ2*sinθ1 + 340*cosθ3*sinθ1*sinθ2 + 73*cosθ5*(cosθ2*sinθ1*sinθ3 + cosθ3*sinθ1*sinθ2) - 45*sinθ1*sinθ2*sinθ3 + 50*sinθ1 + 73*sinθ5*(-cosθ1*sinθ4 + cosθ4*(cosθ2*cosθ3*sinθ1 - sinθ1*sinθ2*sinθ3))], 
	[cosθ6*(cosθ4*cosθ5*(cosθ2*sinθ3 + cosθ3*sinθ2) - sinθ5*(-cosθ2*cosθ3 + sinθ2*sinθ3)) - sinθ4*sinθ6*(cosθ2*sinθ3 + cosθ3*sinθ2), -cosθ6*sinθ4*(cosθ2*sinθ3 + cosθ3*sinθ2) - sinθ6*(cosθ4*cosθ5*(cosθ2*sinθ3 + cosθ3*sinθ2) - sinθ5*(-cosθ2*cosθ3 + sinθ2*sinθ3)), cosθ4*sinθ5*(cosθ2*sinθ3 + cosθ3*sinθ2) + cosθ5*(-cosθ2*cosθ3 + sinθ2*sinθ3), -340*cosθ2*cosθ3 + 45*cosθ2*sinθ3 + 45*cosθ3*sinθ2 + 73*cosθ4*sinθ5*(cosθ2*sinθ3 + cosθ3*sinθ2) + 73*cosθ5*(-cosθ2*cosθ3 + sinθ2*sinθ3) + 340*sinθ2*sinθ3 + 330*sinθ2 + 345], 
	[0, 0, 0, 1]])

	X = 45*cosθ1*cosθ2*cosθ3 + 340*cosθ1*cosθ2*sinθ3 + 330*cosθ1*cosθ2 + 340*cosθ1*cosθ3*sinθ2 - 45*cosθ1*sinθ2*sinθ3 + 50*cosθ1 + 73*cosθ5*(cosθ1*cosθ2*sinθ3 + cosθ1*cosθ3*sinθ2) + 73*sinθ5*(cosθ4*(cosθ1*cosθ2*cosθ3 - cosθ1*sinθ2*sinθ3) + sinθ1*sinθ4)
	Y = 45*cosθ2*cosθ3*sinθ1 + 340*cosθ2*sinθ1*sinθ3 + 330*cosθ2*sinθ1 + 340*cosθ3*sinθ1*sinθ2 + 73*cosθ5*(cosθ2*sinθ1*sinθ3 + cosθ3*sinθ1*sinθ2) - 45*sinθ1*sinθ2*sinθ3 + 50*sinθ1 + 73*sinθ5*(-cosθ1*sinθ4 + cosθ4*(cosθ2*cosθ3*sinθ1 - sinθ1*sinθ2*sinθ3))
	Z = -340*cosθ2*cosθ3 + 45*cosθ2*sinθ3 + 45*cosθ3*sinθ2 + 73*cosθ4*sinθ5*(cosθ2*sinθ3 + cosθ3*sinθ2) + 73*cosθ5*(-cosθ2*cosθ3 + sinθ2*sinθ3) + 340*sinθ2*sinθ3 + 330*sinθ2 + 345

	e_ωθ = R = sympy.Matrix([[r11, r12, r13],
							 [r21, r22, r23],
							 [r31, r32, r33]])

	e_ωθ = sympy.Matrix([ω1*ω1*(1-cosθ)+cosθ,    ω1*ω2*(1-cosθ)-ω3*sinθ, ω1*ω3*(1-cosθ)+ω2*sinθ],
						[ω1*ω2*(1-cosθ)+ω3*sinθ, ω2*ω2*(1-cosθ)+cosθ,    ω2*ω3*(1-cosθ)-ω1*sinθ],
						[ω1*ω3*(1-cosθ)-ω2*sinθ, ω2*ω3*(1-cosθ)+ω1*sinθ, ω3*ω3*(1-cosθ)+cosθ])
	
	θ = arccos((r11 + r22 + r33  - 1)/2)

	RPY(Roll, Pitch, Yaw) = Rz(Roll)*Ry(Pitch)*Rx(Yaw)
                     =[[cos(Roll)*cos(Pitch)  -sin(Roll)*cos(Yaw)+cos(Roll)*sin(Pitch)*sin(Yaw)     sin(Roll)*sin(Yaw)+cos(Roll)*sin(Pitch)*cos(Yaw)]
	 				   [sin(Roll)*cos(Pitch)   cos(Roll)*cos(Yaw)+sin(Roll)*sin(Pitch)*sin(Yaw)    -cos(Roll)*sin(Yaw)+sin(Roll)*sin(Pitch)*cos(Yaw)]
	  				   [-sin(Pitch)               cos(Pitch)*sin(Yaw)                                                    cos(Pitch)*cos(Yaw)        ]]
 
	'''

	'''目标点计算'''
	Roll = Roll*np.pi/180
	Pitch = Pitch*np.pi/180
	Yaw = Yaw*np.pi/180 

	r11 = np.cos(Pitch)*np.cos(Roll)
	r12 = np.cos(Roll)*np.sin(Pitch)*np.sin(Yaw) - np.cos(Yaw)*np.sin(Roll)
	r13 = np.cos(Roll)*np.cos(Yaw)*np.sin(Pitch) + np.sin(Roll)*np.sin(Yaw)

	r21 = np.cos(Pitch)*np.sin(Roll)
	r22 = np.cos(Roll)*np.cos(Yaw) + np.sin(Pitch)*np.sin(Roll)*np.sin(Yaw)
	r23 = -np.cos(Roll)*np.sin(Yaw) + np.cos(Yaw)*np.sin(Pitch)*np.sin(Roll)

	r31 = -np.sin(Pitch)
	r32 = np.cos(Pitch)*np.sin(Yaw)
	r33 = np.cos(Pitch)*np.cos(Yaw)
	'''目标位置'''
	Ts_target = np.array([[r11,r12,r13,X],
						  [r21,r22,r23,Y],
						  [r31,r32,r33,Z],
						  [0,  0,  0,  1]])

	M = np.array([[1, 0, 0, 50+330+45],
				  [0,-1, 0,         0],
				  [0, 0,-1,345-340-73],
				  [0, 0, 0,         1]])

	S1 = np.array([[0],[0],[1],[0],[0],[0]])
	S2 = np.array([[0],[-1],[0],[345],[0],[-50]])
	S3 = np.array([[0],[-1],[0],[345],[0],[-330-50]])
	S4 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])
	S5 = np.array([[0],[-1],[0],[345-340],[0],[-50-330-45]])
	S6 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])

	B1 = np.array([[0],[0],[-1],[0],[-50-330-45],[0]])
	B2 = np.array([[0],[1],[0],[340+73],[0],[-330-45]])
	B3 = np.array([[0],[1],[0],[340+73],[0],[-45]])
	B4 = np.array([[0],[0],[1],[0],[0],[0]])
	B5 = np.array([[0],[1],[0],[73],[0],[0]])	
	B6 = np.array([[0],[0],[1],[0],[0],[0]])

	errorlist = []
	theta1 = []
	theta2 = []
	theta3 = []
	theta4 = []
	theta5 = []
	theta6 = []

	theta = np.array([[theta0_1*np.pi/180],
					  [theta0_2*np.pi/180],
					  [theta0_3*np.pi/180],
					  [theta0_4*np.pi/180],
					  [theta0_5*np.pi/180],					  
					  [theta0_6*np.pi/180]])
	I = np.array([[1,0,0],
				  [0,1,0],
				  [0,0,1]])

	for i in range(steps): 

		e_Sθ1 = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],theta[0][0])
		e_Sθ2 = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],theta[1][0])
		e_Sθ3 = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],theta[2][0])
		e_Sθ4 = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],theta[3][0])
		e_Sθ5 = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],theta[4][0])
		e_Sθ6 = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],theta[5][0])

		Ts_b12 = np.dot(e_Sθ1, e_Sθ2)
		Ts_b13 = np.dot(Ts_b12, e_Sθ3)
		Ts_b14 = np.dot(Ts_b13, e_Sθ4)	
		Ts_b15 = np.dot(Ts_b14, e_Sθ5)
		Ts_b16 = np.dot(Ts_b15, e_Sθ6)
		Ts_b = np.dot(Ts_b16,M) 

		Ts_b_Inverse = np.linalg.inv(Ts_b)		
		Tb_s = Ts_b_Inverse 
		Tb_target = np.dot(Tb_s,Ts_target)

		nx = Tb_target[0,0]
		ny = Tb_target[1,0]
		nz = Tb_target[2,0]
		ox = Tb_target[0,1]
		oy = Tb_target[1,1]
		oz = Tb_target[2,1]
		ax = Tb_target[0,2]
		ay = Tb_target[1,2]
		az = Tb_target[2,2]
		px = Tb_target[0,3]
		py = Tb_target[1,3]
		pz = Tb_target[2,3]

		θi = np. arccos((nx + oy + az  - 1)/2)

		ωi1 = (oz-ay)/(2*np.sin(θi))
		ωi2 = (ax-nz)/(2*np.sin(θi))
		ωi3 = (ny-ox)/(2*np.sin(θi))

		ωi = np.random.random(9)
		ωi = ωi.reshape(3, 3)

		ωi[0,0] = 0               
		ωi[0,1] = -ωi3                   
		ωi[0,2] = ωi2        
		ωi[1,0] = ωi3              
		ωi[1,1] = 0                  
		ωi[1,2] = -ωi1    
		ωi[2,0] = -ωi2               
		ωi[2,1] = ωi1                   
		ωi[2,2] = 0  

		'''
		Ri_Inverse = I/θi-ωi/2+ np.dot((1/θi-1/(np.tan(θi/2)*2)),np.dot(ωi,ωi))

		'''
		Ri = np.dot(I,θi)+np.dot((1-np.cos(θi)),ωi)+np.dot((θi-np.sin(θi)),np.dot(ωi,ωi))			
		Ri_Inverse = np.linalg.inv(Ri)

		Pi = [[px],
			  [py],
			  [pz]]
		vi = np.dot(Ri_Inverse,Pi)

		vi1 = vi[0][0]
		vi2 = vi[1][0]
		vi3 = vi[2][0]

		Vbi = np.array([[ωi1*θi],
						[ωi2*θi],
						[ωi3*θi],
						[vi1*θi],
						[vi2*θi],
						[vi3*θi]])
	
		AdT_b1 = PoE.AdT(B1[0][0],B1[1][0],B1[2][0],B1[3][0],B1[4][0],B1[5][0],-theta[0][0])
		AdT_b2 = PoE.AdT(B2[0][0],B2[1][0],B2[2][0],B2[3][0],B2[4][0],B2[5][0],-theta[1][0])
		AdT_b3 = PoE.AdT(B3[0][0],B3[1][0],B3[2][0],B3[3][0],B3[4][0],B3[5][0],-theta[2][0])
		AdT_b4 = PoE.AdT(B4[0][0],B4[1][0],B4[2][0],B4[3][0],B4[4][0],B4[5][0],-theta[3][0])
		AdT_b5 = PoE.AdT(B5[0][0],B5[1][0],B5[2][0],B5[3][0],B5[4][0],B5[5][0],-theta[4][0])
		AdT_b6 = PoE.AdT(B6[0][0],B6[1][0],B6[2][0],B6[3][0],B6[4][0],B6[5][0],-theta[5][0])

		AdT_b65 = np.dot(AdT_b6,AdT_b5)		
		AdT_b654= np.dot(AdT_b65,AdT_b4)
		AdT_b6543 = np.dot(AdT_b654,AdT_b3)
		AdT_b65432 = np.dot(AdT_b6543,AdT_b2)

		Jb1 = np.dot(AdT_b65432,B1)
		Jb2 = np.dot(AdT_b6543,B2)
		Jb3 = np.dot(AdT_b654,B3)
		Jb4 = np.dot(AdT_b65,B4)
		Jb5 = np.dot(AdT_b6,B5)
		Jb6 = B6

		Jb = np.random.random(36)
		Jb = Jb.reshape(6, 6)

		Jb[0,0] = Jb1[0][0]
		Jb[1,0] = Jb1[1][0]
		Jb[2,0] = Jb1[2][0]
		Jb[3,0] = Jb1[3][0]
		Jb[4,0] = Jb1[4][0]
		Jb[5,0] = Jb1[5][0]

		Jb[0,1] = Jb2[0][0]
		Jb[1,1] = Jb2[1][0]
		Jb[2,1] = Jb2[2][0]
		Jb[3,1] = Jb2[3][0]
		Jb[4,1] = Jb2[4][0]
		Jb[5,1] = Jb2[5][0]

		Jb[0,2] = Jb3[0][0]
		Jb[1,2] = Jb3[1][0]
		Jb[2,2] = Jb3[2][0]
		Jb[3,2] = Jb3[3][0]
		Jb[4,2] = Jb3[4][0]
		Jb[5,2] = Jb3[5][0]

		Jb[0,3] = Jb4[0][0]
		Jb[1,3] = Jb4[1][0]
		Jb[2,3] = Jb4[2][0]
		Jb[3,3] = Jb4[3][0]
		Jb[4,3] = Jb4[4][0]
		Jb[5,3] = Jb4[5][0]

		Jb[0,4] = Jb5[0][0]
		Jb[1,4] = Jb5[1][0]
		Jb[2,4] = Jb5[2][0]
		Jb[3,4] = Jb5[3][0]
		Jb[4,4] = Jb5[4][0]
		Jb[5,4] = Jb5[5][0]

		Jb[0,5] = Jb6[0][0]
		Jb[1,5] = Jb6[1][0]
		Jb[2,5] = Jb6[2][0]
		Jb[3,5] = Jb6[3][0]
		Jb[4,5] = Jb6[4][0]
		Jb[5,5] = Jb6[5][0]

		'''求矩阵的逆'''
		Jb_Inverse = np.linalg.inv(Jb)

		thetaX = theta
		theta = theta + np.dot(Jb_Inverse,Vbi)
		deltatheta = theta-thetaX

		#errorlist.append(np.linalg.norm(deltatheta))
		e = np.sqrt(np.dot(deltatheta[0][0],deltatheta[0][0]) + np.dot(deltatheta[1][0],deltatheta[1][0]) + np.dot(deltatheta[2][0],deltatheta[2][0]) + np.dot(deltatheta[3][0],deltatheta[3][0])+ np.dot(deltatheta[4][0],deltatheta[4][0])+ np.dot(deltatheta[5][0],deltatheta[5][0]))
		errorlist.append(e)
			
		theta1.append(theta[0][0]*180/np.pi)
		theta2.append(theta[1][0]*180/np.pi)
		theta3.append(theta[2][0]*180/np.pi)
		theta4.append(theta[3][0]*180/np.pi)
		theta5.append(theta[4][0]*180/np.pi)
		theta6.append(theta[5][0]*180/np.pi)

		if errorlist[-1] < error:
			print('iteration: ', theta*180/np.pi%360)
			
			break
	x_axis = [i for i in range(len(errorlist))]
	print('迭代次数 = ',i+1)
	plt.figure()
	plt.plot(x_axis, errorlist)
	plt.plot(x_axis, theta1)	
	plt.plot(x_axis, theta2)
	plt.plot(x_axis, theta3)
	plt.plot(x_axis, theta4)
	plt.plot(x_axis, theta5)
	plt.plot(x_axis, theta6)
	plt.show()

	return [theta[0][0]*180/np.pi%360,theta[1][0]*180/np.pi%360,theta[2][0]*180/np.pi%360,theta[3][0]*180/np.pi%360,theta[4][0]*180/np.pi%360,theta[5][0]*180/np.pi%360]

def MZ07_Trajectory_Planning_StraightLine_TriplePolynomial(θ1_origin, θ2_origin,θ3_origin,θ4_origin,θ5_origin,θ6_origin,θ1_target, θ2_target,θ3_target,θ4_target,θ5_target,θ6_target,T=3,t=0.01):


	M = np.array([[1, 0, 0, 50+330+45],
				  [0,-1, 0,         0],
			 	  [0, 0,-1,345-340-73],
				  [0, 0, 0,         1]])

	S1 = np.array([[0],[0],[1],[0],[0],[0]])
	S2 = np.array([[0],[-1],[0],[345],[0],[-50]])
	S3 = np.array([[0],[-1],[0],[345],[0],[-330-50]])
	S4 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])
	S5 = np.array([[0],[-1],[0],[345-340],[0],[-50-330-45]])
	S6 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])
	
	'''设定终点'''
	e_Sθ1_target = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ1_target*np.pi/180)
	e_Sθ2_target = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ2_target*np.pi/180)
	e_Sθ3_target = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ3_target*np.pi/180)
	e_Sθ4_target = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ4_target*np.pi/180)
	e_Sθ5_target = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],θ5_target*np.pi/180)
	e_Sθ6_target = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],θ6_target*np.pi/180)

	T12_target = np.dot(e_Sθ1_target, e_Sθ2_target)
	T13_target = np.dot(T12_target, e_Sθ3_target)
	T14_target = np.dot(T13_target, e_Sθ4_target)
	T15_target = np.dot(T14_target, e_Sθ5_target)
	T16_target = np.dot(T15_target, e_Sθ6_target)
	Ts_target = np.dot(T16_target, M) 

	'''设定起点点'''
	e_Sθ1_origin = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ1_origin*np.pi/180)
	e_Sθ2_origin = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ2_origin*np.pi/180)
	e_Sθ3_origin = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ3_origin*np.pi/180)
	e_Sθ4_origin = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ4_origin*np.pi/180)
	e_Sθ5_origin = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],θ5_origin*np.pi/180)
	e_Sθ6_origin = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],θ6_origin*np.pi/180)

	T12_origin = np.dot(e_Sθ1_origin, e_Sθ2_origin)
	T13_origin = np.dot(T12_origin, e_Sθ3_origin)
	T14_origin = np.dot(T13_origin, e_Sθ4_origin)
	T15_origin = np.dot(T14_origin, e_Sθ5_origin)
	T16_origin = np.dot(T15_origin, e_Sθ6_origin)
	Ts_origin = np.dot(T16_origin, M) 

	px = []
	py = []
	pz = []
	vx = []
	vy = []
	vz = []
	v = []
	a0 = 0 
	a1 = 0
	a2 = 3/(T*T)
	a3 = -2/(T*T*T)
	for i in np.arange(0,T+t,t):
		s = a0 + a1*i+a2*i*i+a3*i*i*i
		Tst = Ts_origin + s*(Ts_target-Ts_origin)

		px.append(Tst[0,3])
		py.append(Tst[1,3])
		pz.append(Tst[2,3])

		V_Tst = (6*i/(T*T)-6*i*i/(T*T*T))*(Ts_target-Ts_origin)

		vx.append(V_Tst[0,3])
		vy.append(V_Tst[1,3])
		vz.append(V_Tst[2,3])

		v.append(np.sqrt(vx[-1]*vx[-1] + vy[-1]*vy[-1] + vz[-1]*vz[-1]))
		
	print(v)
	x_axis = [i for i in range(len(px))]
	plt.figure()
	#plta = plt.gca(projection='3d')
	#plta.plot(px, py, pz)
	plt.plot(x_axis,v)
	plt.show()
def MZ07_Trajectory_Planning_StraightLine_QuinticPolynomial(θ1_origin, θ2_origin,θ3_origin,θ4_origin,θ5_origin,θ6_origin,θ1_target, θ2_target,θ3_target,θ4_target,θ5_target,θ6_target,T=3,t=0.01):
	a0 = 0 
	a1 = 0
	a2 = 0
	a3 = 10/(T*T*T)
	a4 = -15/(T*T*T*T)
	a5 = 6/(T*T*T*T*T)

	M = np.array([[1, 0, 0, 50+330+45],
				  [0,-1, 0,         0],
			 	  [0, 0,-1,345-340-73],
				  [0, 0, 0,         1]])

	S1 = np.array([[0],[0],[1],[0],[0],[0]])
	S2 = np.array([[0],[-1],[0],[345],[0],[-50]])
	S3 = np.array([[0],[-1],[0],[345],[0],[-330-50]])
	S4 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])
	S5 = np.array([[0],[-1],[0],[345-340],[0],[-50-330-45]])
	S6 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])
	
	'''设定终点'''
	e_Sθ1_target = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ1_target*np.pi/180)
	e_Sθ2_target = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ2_target*np.pi/180)
	e_Sθ3_target = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ3_target*np.pi/180)
	e_Sθ4_target = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ4_target*np.pi/180)
	e_Sθ5_target = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],θ5_target*np.pi/180)
	e_Sθ6_target = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],θ6_target*np.pi/180)

	T12_target = np.dot(e_Sθ1_target, e_Sθ2_target)
	T13_target = np.dot(T12_target, e_Sθ3_target)
	T14_target = np.dot(T13_target, e_Sθ4_target)
	T15_target = np.dot(T14_target, e_Sθ5_target)
	T16_target = np.dot(T15_target, e_Sθ6_target)
	Ts_target = np.dot(T16_target, M) 

	'''设定起点点'''
	e_Sθ1_origin = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ1_origin*np.pi/180)
	e_Sθ2_origin = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ2_origin*np.pi/180)
	e_Sθ3_origin = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ3_origin*np.pi/180)
	e_Sθ4_origin = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ4_origin*np.pi/180)
	e_Sθ5_origin = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],θ5_origin*np.pi/180)
	e_Sθ6_origin = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],θ6_origin*np.pi/180)

	T12_origin = np.dot(e_Sθ1_origin, e_Sθ2_origin)
	T13_origin = np.dot(T12_origin, e_Sθ3_origin)
	T14_origin = np.dot(T13_origin, e_Sθ4_origin)
	T15_origin = np.dot(T14_origin, e_Sθ5_origin)
	T16_origin = np.dot(T15_origin, e_Sθ6_origin)
	Ts_origin = np.dot(T16_origin, M) 

	px = []
	py = []
	pz = []
	vx = []
	vy = []
	vz = []
	v = []
	
	for i in np.arange(0,T+t,t):
		s = a0+a1*i+a2*i*i+a3*i*i*i+a4*i*i*i*i+a5*i*i*i*i*i
		Tst = Ts_origin + s*(Ts_target-Ts_origin)

		px.append(Tst[0,3])
		py.append(Tst[1,3])
		pz.append(Tst[2,3])

		V_Tst = (6*i/(T*T)-6*i*i/(T*T*T))*(Ts_target-Ts_origin)

		vx.append(V_Tst[0,3])
		vy.append(V_Tst[1,3])
		vz.append(V_Tst[2,3])

		v.append(np.sqrt(vx[-1]*vx[-1] + vy[-1]*vy[-1] + vz[-1]*vz[-1]))
		
	print(v)
	x_axis = [i for i in range(len(px))]
	plt.figure()
	#plta = plt.gca(projection='3d')
	#plta.plot(px, py, pz)
	plt.plot(x_axis,v)
	plt.show()
def MZ07_Trajectory_Planning_StraightLine_TrapezoidalMotion_Parameter_v_T(θ1_origin, θ2_origin,θ3_origin,θ4_origin,θ5_origin,θ6_origin,θ1_target, θ2_target,θ3_target,θ4_target,θ5_target,θ6_target,v = 10, T=0.15,t=0.01):
	'''最大速度v，加速度a(不控制),时间间隔t,运动时间T'''
	'''梯形图的必要条件2>=vT>1'''
	M = np.array([[1, 0, 0, 50+330+45],
				  [0,-1, 0,         0],
			 	  [0, 0,-1,345-340-73],
				  [0, 0, 0,         1]])

	S1 = np.array([[0],[0],[1],[0],[0],[0]])
	S2 = np.array([[0],[-1],[0],[345],[0],[-50]])
	S3 = np.array([[0],[-1],[0],[345],[0],[-330-50]])
	S4 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])
	S5 = np.array([[0],[-1],[0],[345-340],[0],[-50-330-45]])
	S6 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])
	
	'''设定终点'''
	e_Sθ1_target = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ1_target*np.pi/180)
	e_Sθ2_target = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ2_target*np.pi/180)
	e_Sθ3_target = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ3_target*np.pi/180)
	e_Sθ4_target = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ4_target*np.pi/180)
	e_Sθ5_target = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],θ5_target*np.pi/180)
	e_Sθ6_target = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],θ6_target*np.pi/180)

	T12_target = np.dot(e_Sθ1_target, e_Sθ2_target)
	T13_target = np.dot(T12_target, e_Sθ3_target)
	T14_target = np.dot(T13_target, e_Sθ4_target)
	T15_target = np.dot(T14_target, e_Sθ5_target)
	T16_target = np.dot(T15_target, e_Sθ6_target)
	Ts_target = np.dot(T16_target, M) 

	'''设定起点'''
	e_Sθ1_origin = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ1_origin*np.pi/180)
	e_Sθ2_origin = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ2_origin*np.pi/180)
	e_Sθ3_origin = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ3_origin*np.pi/180)
	e_Sθ4_origin = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ4_origin*np.pi/180)
	e_Sθ5_origin = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],θ5_origin*np.pi/180)
	e_Sθ6_origin = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],θ6_origin*np.pi/180)

	T12_origin = np.dot(e_Sθ1_origin, e_Sθ2_origin)
	T13_origin = np.dot(T12_origin, e_Sθ3_origin)
	T14_origin = np.dot(T13_origin, e_Sθ4_origin)
	T15_origin = np.dot(T14_origin, e_Sθ5_origin)
	T16_origin = np.dot(T15_origin, e_Sθ6_origin)
	Ts_origin = np.dot(T16_origin, M) 

	px = []
	py = []
	pz = []
	vx = []
	vy = []
	vz = []
	V = []

	a = v*v/(v*T-1)

	for i in np.arange(0,T+t,t):
		
		if i>=0 and i <=v/a:
			
			Tst = Ts_origin + (a*i*i/2)*(Ts_target-Ts_origin)
			V_Tst = (a*i)*(Ts_target-Ts_origin)

		if i>v/a and i <=T-v/a:
			Tst = Ts_origin + (v*i-v*v/(2*a))*(Ts_target-Ts_origin)
			V_Tst = (v)*(Ts_target-Ts_origin)			

		if i>T-v/a and i <=T:
			Tst = Ts_origin + ((2*a*v*T-2*v*v-a*a*(i-T)*(i-T))/(2*a))*(Ts_target-Ts_origin)
			V_Tst = (a*(T-i))*(Ts_target-Ts_origin)	



		px.append(Tst[0,3])
		py.append(Tst[1,3])
		pz.append(Tst[2,3])

		
		vx.append(V_Tst[0,3])
		vy.append(V_Tst[1,3])
		vz.append(V_Tst[2,3])

		V.append(np.sqrt(vx[-1]*vx[-1] + vy[-1]*vy[-1] + vz[-1]*vz[-1]))
		print(i)
		
	print(V)
	x_axis = [i for i in range(len(px))]
	plt.figure()
	plta = plt.gca(projection='3d')
	plta.plot(px, py, pz)
	#plt.plot(x_axis,V)
	plta.show()
def MZ07_Trajectory_Planning_StraightLine_TrapezoidalMotion_Parameter_v_a(θ1_origin, θ2_origin,θ3_origin,θ4_origin,θ5_origin,θ6_origin,θ1_target, θ2_target,θ3_target,θ4_target,θ5_target,θ6_target,v = 10, a=150,t=0.0001):
	'''最大速度v，加速度a,时间间隔t,运动时间T(不控制)'''
	'''梯形图的必要条件v*v/a<=1'''
	M = np.array([[1, 0, 0, 50+330+45],
				  [0,-1, 0,         0],
			 	  [0, 0,-1,345-340-73],
				  [0, 0, 0,         1]])

	S1 = np.array([[0],[0],[1],[0],[0],[0]])
	S2 = np.array([[0],[-1],[0],[345],[0],[-50]])
	S3 = np.array([[0],[-1],[0],[345],[0],[-330-50]])
	S4 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])
	S5 = np.array([[0],[-1],[0],[345-340],[0],[-50-330-45]])
	S6 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])
	
	'''设定终点'''
	e_Sθ1_target = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ1_target*np.pi/180)
	e_Sθ2_target = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ2_target*np.pi/180)
	e_Sθ3_target = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ3_target*np.pi/180)
	e_Sθ4_target = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ4_target*np.pi/180)
	e_Sθ5_target = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],θ5_target*np.pi/180)
	e_Sθ6_target = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],θ6_target*np.pi/180)

	T12_target = np.dot(e_Sθ1_target, e_Sθ2_target)
	T13_target = np.dot(T12_target, e_Sθ3_target)
	T14_target = np.dot(T13_target, e_Sθ4_target)
	T15_target = np.dot(T14_target, e_Sθ5_target)
	T16_target = np.dot(T15_target, e_Sθ6_target)
	Ts_target = np.dot(T16_target, M) 

	'''设定起点'''
	e_Sθ1_origin = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ1_origin*np.pi/180)
	e_Sθ2_origin = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ2_origin*np.pi/180)
	e_Sθ3_origin = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ3_origin*np.pi/180)
	e_Sθ4_origin = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ4_origin*np.pi/180)
	e_Sθ5_origin = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],θ5_origin*np.pi/180)
	e_Sθ6_origin = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],θ6_origin*np.pi/180)

	T12_origin = np.dot(e_Sθ1_origin, e_Sθ2_origin)
	T13_origin = np.dot(T12_origin, e_Sθ3_origin)
	T14_origin = np.dot(T13_origin, e_Sθ4_origin)
	T15_origin = np.dot(T14_origin, e_Sθ5_origin)
	T16_origin = np.dot(T15_origin, e_Sθ6_origin)
	Ts_origin = np.dot(T16_origin, M) 

	px = []
	py = []
	pz = []
	vx = []
	vy = []
	vz = []
	V = []

	T = (a+v*v)/(v*a)

	for i in np.arange(0,T+t,t):
		
		if i>=0 and i <=v/a:
			
			Tst = Ts_origin + (a*i*i/2)*(Ts_target-Ts_origin)
			V_Tst = (a*i)*(Ts_target-Ts_origin)

		if i>v/a and i <=T-v/a:
			Tst = Ts_origin + (v*i-v*v/(2*a))*(Ts_target-Ts_origin)
			V_Tst = (v)*(Ts_target-Ts_origin)			

		if i>T-v/a and i <=T:
			Tst = Ts_origin + ((2*a*v*T-2*v*v-a*a*(i-T)*(i-T))/(2*a))*(Ts_target-Ts_origin)
			V_Tst = (a*(T-i))*(Ts_target-Ts_origin)	

		px.append(Tst[0,3])
		py.append(Tst[1,3])
		pz.append(Tst[2,3])

		
		vx.append(V_Tst[0,3])
		vy.append(V_Tst[1,3])
		vz.append(V_Tst[2,3])

		V.append(np.sqrt(vx[-1]*vx[-1] + vy[-1]*vy[-1] + vz[-1]*vz[-1]))
		print(i)
		
	print(V)
	x_axis = [i for i in range(len(px))]
	plt.figure()
	#plta = plt.gca(projection='3d')
	#plta.plot(px, py, pz)
	plt.plot(x_axis,V)
	plt.show()
def MZ07_Trajectory_Planning_StraightLine_TrapezoidalMotion_Parameter_a_T(θ1_origin, θ2_origin,θ3_origin,θ4_origin,θ5_origin,θ6_origin,θ1_target, θ2_target,θ3_target,θ4_target,θ5_target,θ6_target,a = 10, T =1,t=0.0001):
	'''最大速度v(不控制)，加速度a,时间间隔t,运动时间T'''
	'''梯形图的必要条件a*T*T>=4'''
	M = np.array([[1, 0, 0, 50+330+45],
				  [0,-1, 0,         0],
			 	  [0, 0,-1,345-340-73],
				  [0, 0, 0,         1]])

	S1 = np.array([[0],[0],[1],[0],[0],[0]])
	S2 = np.array([[0],[-1],[0],[345],[0],[-50]])
	S3 = np.array([[0],[-1],[0],[345],[0],[-330-50]])
	S4 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])
	S5 = np.array([[0],[-1],[0],[345-340],[0],[-50-330-45]])
	S6 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])
	
	'''设定终点'''
	e_Sθ1_target = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ1_target*np.pi/180)
	e_Sθ2_target = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ2_target*np.pi/180)
	e_Sθ3_target = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ3_target*np.pi/180)
	e_Sθ4_target = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ4_target*np.pi/180)
	e_Sθ5_target = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],θ5_target*np.pi/180)
	e_Sθ6_target = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],θ6_target*np.pi/180)

	T12_target = np.dot(e_Sθ1_target, e_Sθ2_target)
	T13_target = np.dot(T12_target, e_Sθ3_target)
	T14_target = np.dot(T13_target, e_Sθ4_target)
	T15_target = np.dot(T14_target, e_Sθ5_target)
	T16_target = np.dot(T15_target, e_Sθ6_target)
	Ts_target = np.dot(T16_target, M) 

	'''设定起点'''
	e_Sθ1_origin = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ1_origin*np.pi/180)
	e_Sθ2_origin = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ2_origin*np.pi/180)
	e_Sθ3_origin = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ3_origin*np.pi/180)
	e_Sθ4_origin = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ4_origin*np.pi/180)
	e_Sθ5_origin = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],θ5_origin*np.pi/180)
	e_Sθ6_origin = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],θ6_origin*np.pi/180)

	T12_origin = np.dot(e_Sθ1_origin, e_Sθ2_origin)
	T13_origin = np.dot(T12_origin, e_Sθ3_origin)
	T14_origin = np.dot(T13_origin, e_Sθ4_origin)
	T15_origin = np.dot(T14_origin, e_Sθ5_origin)
	T16_origin = np.dot(T15_origin, e_Sθ6_origin)
	Ts_origin = np.dot(T16_origin, M) 

	px = []
	py = []
	pz = []
	vx = []
	vy = []
	vz = []
	V = []

	v = (a*T-np.sqrt(a)*np.sqrt(a*T*T-4))/2

	for i in np.arange(0,T+t,t):
		
		if i>=0 and i <=v/a:
			
			Tst = Ts_origin + (a*i*i/2)*(Ts_target-Ts_origin)
			V_Tst = (a*i)*(Ts_target-Ts_origin)

		if i>v/a and i <=T-v/a:
			Tst = Ts_origin + (v*i-v*v/(2*a))*(Ts_target-Ts_origin)
			V_Tst = (v)*(Ts_target-Ts_origin)			

		if i>T-v/a and i <=T:
			Tst = Ts_origin + ((2*a*v*T-2*v*v-a*a*(i-T)*(i-T))/(2*a))*(Ts_target-Ts_origin)
			V_Tst = (a*(T-i))*(Ts_target-Ts_origin)	

		px.append(Tst[0,3])
		py.append(Tst[1,3])
		pz.append(Tst[2,3])

		
		vx.append(V_Tst[0,3])
		vy.append(V_Tst[1,3])
		vz.append(V_Tst[2,3])

		V.append(np.sqrt(vx[-1]*vx[-1] + vy[-1]*vy[-1] + vz[-1]*vz[-1]))
		print(i)
		
	print(V)
	x_axis = [i for i in range(len(px))]
	plt.figure()
	plta = plt.gca(projection='3d')
	#plta.plot(px, py, pz)
	plt.plot(x_axis,V)
	plt.show()

def MZ07_Trajectory_Planning_StraightLine_TriplePolynomial_AngleOutput(θ1_origin, θ2_origin,θ3_origin,θ4_origin,θ5_origin,θ6_origin,θ1_target, θ2_target,θ3_target,θ4_target,θ5_target,θ6_target,T=3,t=0.01,steps=10,error=0.001):
	a0 = 0 
	a1 = 0
	a2 = 3/(T*T)
	a3 = -2/(T*T*T)

	M = np.array([[1, 0, 0, 50+330+45],
				  [0,-1, 0,         0],
			 	  [0, 0,-1,345-340-73],
				  [0, 0, 0,         1]])
	
	errorlist = []
	theta1 = []
	theta2 = []
	theta3 = []
	theta4 = []
	theta5 = []
	theta6 = []

	I = np.array([[1,0,0],
				  [0,1,0],
				  [0,0,1]])

	S1 = np.array([[0],[0],[1],[0],[0],[0]])
	S2 = np.array([[0],[-1],[0],[345],[0],[-50]])
	S3 = np.array([[0],[-1],[0],[345],[0],[-330-50]])
	S4 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])
	S5 = np.array([[0],[-1],[0],[345-340],[0],[-50-330-45]])
	S6 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])

	B1 = np.array([[0],[0],[-1],[0],[-50-330-45],[0]])
	B2 = np.array([[0],[1],[0],[340+73],[0],[-330-45]])
	B3 = np.array([[0],[1],[0],[340+73],[0],[-45]])
	B4 = np.array([[0],[0],[1],[0],[0],[0]])
	B5 = np.array([[0],[1],[0],[73],[0],[0]])	
	B6 = np.array([[0],[0],[1],[0],[0],[0]])

	'''设定终点'''
	θ_target = np.array([[θ1_target*np.pi/180],
					  	 [θ2_target*np.pi/180],
						 [θ3_target*np.pi/180],
						 [θ4_target*np.pi/180],
						 [θ5_target*np.pi/180],
						 [θ6_target*np.pi/180]])

	e_Sθ1_target = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ_target[0][0])
	e_Sθ2_target = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ_target[1][0])
	e_Sθ3_target = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ_target[2][0])
	e_Sθ4_target = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ_target[3][0])
	e_Sθ5_target = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],θ_target[4][0])
	e_Sθ6_target = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],θ_target[5][0])

	T12_target = np.dot(e_Sθ1_target, e_Sθ2_target)
	T13_target = np.dot(T12_target, e_Sθ3_target)
	T14_target = np.dot(T13_target, e_Sθ4_target)
	T15_target = np.dot(T14_target, e_Sθ5_target)
	T16_target = np.dot(T15_target, e_Sθ6_target)
	Ts_target = np.dot(T16_target, M) 

	'''设定起点'''
	θ_origin = np.array([[θ1_origin*np.pi/180],
						 [θ2_origin*np.pi/180],
						 [θ3_origin*np.pi/180],
						 [θ4_origin*np.pi/180],
					 	 [θ5_origin*np.pi/180],
						 [θ6_origin*np.pi/180]])

	e_Sθ1_origin = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ_origin[0][0])
	e_Sθ2_origin = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ_origin[1][0])
	e_Sθ3_origin = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ_origin[2][0])
	e_Sθ4_origin = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ_origin[3][0])
	e_Sθ5_origin = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],θ_origin[4][0])
	e_Sθ6_origin = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],θ_origin[5][0])

	T12_origin = np.dot(e_Sθ1_origin, e_Sθ2_origin)
	T13_origin = np.dot(T12_origin, e_Sθ3_origin)
	T14_origin = np.dot(T13_origin, e_Sθ4_origin)
	T15_origin = np.dot(T14_origin, e_Sθ5_origin)
	T16_origin = np.dot(T15_origin, e_Sθ6_origin)
	Ts_origin = np.dot(T16_origin, M) 

	Px_Path = []
	Py_Path = []
	Pz_Path = []
	vx_Path = []
	vy_Path = []
	vz_Path = []
	v_Path = []
	for i in np.arange(0,T+t,t):
		if i == 0:

			theta = θ_origin
			theta1.append(theta[0][0]*180/np.pi)
			theta2.append(theta[1][0]*180/np.pi)
			theta3.append(theta[2][0]*180/np.pi)
			theta4.append(theta[3][0]*180/np.pi)
			theta5.append(theta[4][0]*180/np.pi)
			theta6.append(theta[5][0]*180/np.pi)

			Px_Path.append(Ts_origin[0,3])
			Py_Path.append(Ts_origin[1,3])
			Pz_Path.append(Ts_origin[2,3])

			V_Tst = (6*i/(T*T)-6*i*i/(T*T*T))*(Ts_target-Ts_origin)

			vx_Path.append(V_Tst[0,3])
			vy_Path.append(V_Tst[1,3])
			vz_Path.append(V_Tst[2,3])

			v_Path.append(np.sqrt(vx_Path[-1]*vx_Path[-1] + vy_Path[-1]*vy_Path[-1] + vz_Path[-1]*vz_Path[-1]))
			
			i = i + t

		Tst = Ts_origin + (3*i*i/(T*T)-2*i*i*i/(T*T*T))*(Ts_target-Ts_origin)

		'''逆运动学数值解'''
		for j in range(steps): 

			e_Sθ1 = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],theta[0][0])
			e_Sθ2 = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],theta[1][0])
			e_Sθ3 = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],theta[2][0])
			e_Sθ4 = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],theta[3][0])
			e_Sθ5 = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],theta[4][0])
			e_Sθ6 = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],theta[5][0])

			Ts_b12 = np.dot(e_Sθ1, e_Sθ2)
			Ts_b13 = np.dot(Ts_b12, e_Sθ3)
			Ts_b14 = np.dot(Ts_b13, e_Sθ4)	
			Ts_b15 = np.dot(Ts_b14, e_Sθ5)
			Ts_b16 = np.dot(Ts_b15, e_Sθ6)
			Ts_b = np.dot(Ts_b16,M) 
			#print('Ts_b =', Ts_b)
			Ts_b_Inverse = np.linalg.inv(Ts_b)		
			Tb_s = Ts_b_Inverse 
			Tb_target = np.dot(Tb_s,Tst)
			
			#print('Tb_target =', Tb_target)

			nx = Tb_target[0,0]
			ny = Tb_target[1,0]
			nz = Tb_target[2,0]
			ox = Tb_target[0,1]
			oy = Tb_target[1,1]
			oz = Tb_target[2,1]
			ax = Tb_target[0,2]
			ay = Tb_target[1,2]
			az = Tb_target[2,2]
			px = Tb_target[0,3]
			py = Tb_target[1,3]
			pz = Tb_target[2,3]

			θi = np. arccos((nx + oy + az  - 1)/2)


			ωi1 = (oz-ay)/(2*np.sin(θi))
			ωi2 = (ax-nz)/(2*np.sin(θi))
			ωi3 = (ny-ox)/(2*np.sin(θi))

			ωi = np.random.random(9)
			ωi = ωi.reshape(3, 3)

			ωi[0,0] = 0               
			ωi[0,1] = -ωi3                   
			ωi[0,2] = ωi2        
			ωi[1,0] = ωi3              
			ωi[1,1] = 0                  
			ωi[1,2] = -ωi1    
			ωi[2,0] = -ωi2               
			ωi[2,1] = ωi1                   
			ωi[2,2] = 0  

			#print('ωi =',ωi)
			'''
			Ri_Inverse = I/θi-ωi/2+ np.dot((1/θi-1/(np.tan(θi/2)*2)),np.dot(ωi,ωi))

			'''
			Ri = np.dot(I,θi)+np.dot((1-np.cos(θi)),ωi)+np.dot((θi-np.sin(θi)),np.dot(ωi,ωi))			

			Ri_Inverse = np.linalg.inv(Ri)

			Pi = [[px],
				  [py],
				  [pz]]
			vi = np.dot(Ri_Inverse,Pi)

			vi1 = vi[0][0]
			vi2 = vi[1][0]
			vi3 = vi[2][0]

			Vbi = np.array([[ωi1*θi],
							[ωi2*θi],
							[ωi3*θi],
							[vi1*θi],
							[vi2*θi],
							[vi3*θi]])
		
			AdT_b1 = PoE.AdT(B1[0][0],B1[1][0],B1[2][0],B1[3][0],B1[4][0],B1[5][0],-theta[0][0])
			AdT_b2 = PoE.AdT(B2[0][0],B2[1][0],B2[2][0],B2[3][0],B2[4][0],B2[5][0],-theta[1][0])
			AdT_b3 = PoE.AdT(B3[0][0],B3[1][0],B3[2][0],B3[3][0],B3[4][0],B3[5][0],-theta[2][0])
			AdT_b4 = PoE.AdT(B4[0][0],B4[1][0],B4[2][0],B4[3][0],B4[4][0],B4[5][0],-theta[3][0])
			AdT_b5 = PoE.AdT(B5[0][0],B5[1][0],B5[2][0],B5[3][0],B5[4][0],B5[5][0],-theta[4][0])
			AdT_b6 = PoE.AdT(B6[0][0],B6[1][0],B6[2][0],B6[3][0],B6[4][0],B6[5][0],-theta[5][0])


			AdT_b65 = np.dot(AdT_b6,AdT_b5)		
			AdT_b654= np.dot(AdT_b65,AdT_b4)
			AdT_b6543 = np.dot(AdT_b654,AdT_b3)
			AdT_b65432 = np.dot(AdT_b6543,AdT_b2)


			Jb1 = np.dot(AdT_b65432,B1)
			Jb2 = np.dot(AdT_b6543,B2)
			Jb3 = np.dot(AdT_b654,B3)
			Jb4 = np.dot(AdT_b65,B4)
			Jb5 = np.dot(AdT_b6,B5)
			Jb6 = B6

			Jb = np.random.random(36)
			Jb = Jb.reshape(6, 6)

			Jb[0,0] = Jb1[0][0]
			Jb[1,0] = Jb1[1][0]
			Jb[2,0] = Jb1[2][0]
			Jb[3,0] = Jb1[3][0]
			Jb[4,0] = Jb1[4][0]
			Jb[5,0] = Jb1[5][0]

			Jb[0,1] = Jb2[0][0]
			Jb[1,1] = Jb2[1][0]
			Jb[2,1] = Jb2[2][0]
			Jb[3,1] = Jb2[3][0]
			Jb[4,1] = Jb2[4][0]
			Jb[5,1] = Jb2[5][0]

			Jb[0,2] = Jb3[0][0]
			Jb[1,2] = Jb3[1][0]
			Jb[2,2] = Jb3[2][0]
			Jb[3,2] = Jb3[3][0]
			Jb[4,2] = Jb3[4][0]
			Jb[5,2] = Jb3[5][0]

			Jb[0,3] = Jb4[0][0]
			Jb[1,3] = Jb4[1][0]
			Jb[2,3] = Jb4[2][0]
			Jb[3,3] = Jb4[3][0]
			Jb[4,3] = Jb4[4][0]
			Jb[5,3] = Jb4[5][0]

			Jb[0,4] = Jb5[0][0]
			Jb[1,4] = Jb5[1][0]
			Jb[2,4] = Jb5[2][0]
			Jb[3,4] = Jb5[3][0]
			Jb[4,4] = Jb5[4][0]
			Jb[5,4] = Jb5[5][0]

			Jb[0,5] = Jb6[0][0]
			Jb[1,5] = Jb6[1][0]
			Jb[2,5] = Jb6[2][0]
			Jb[3,5] = Jb6[3][0]
			Jb[4,5] = Jb6[4][0]
			Jb[5,5] = Jb6[5][0]

			'''求矩阵的逆'''
			Jb_Inverse = np.linalg.inv(Jb)

			thetaX = theta
			theta = theta + np.dot(Jb_Inverse,Vbi)
			deltatheta = theta-thetaX

			#errorlist.append(np.linalg.norm(deltatheta))
			e = np.sqrt(np.dot(deltatheta[0][0],deltatheta[0][0]) + np.dot(deltatheta[1][0],deltatheta[1][0]) + np.dot(deltatheta[2][0],deltatheta[2][0]) + np.dot(deltatheta[3][0],deltatheta[3][0])+ np.dot(deltatheta[4][0],deltatheta[4][0])+ np.dot(deltatheta[5][0],deltatheta[5][0]))
			errorlist.append(e)
				
			if errorlist[-1] < error:
				# print('iteration: ', theta*180/np.pi%360)
				print('迭代次数:',j+1)
				theta1.append(theta[0][0]*180/np.pi)
				theta2.append(theta[1][0]*180/np.pi)
				theta3.append(theta[2][0]*180/np.pi)
				theta4.append(theta[3][0]*180/np.pi)
				theta5.append(theta[4][0]*180/np.pi)
				theta6.append(theta[5][0]*180/np.pi)
				break


		Px_Path.append(Tst[0,3])
		Py_Path.append(Tst[1,3])
		Pz_Path.append(Tst[2,3])


		V_Tst = (6*i/(T*T)-6*i*i/(T*T*T))*(Ts_target-Ts_origin)

		vx_Path.append(V_Tst[0,3])
		vy_Path.append(V_Tst[1,3])
		vz_Path.append(V_Tst[2,3])

		v_Path.append(np.sqrt(vx_Path[-1]*vx_Path[-1] + vy_Path[-1]*vy_Path[-1] + vz_Path[-1]*vz_Path[-1]))


	
	x_axis = [i for i in range(len(theta1))]
	x1_axis = [i for i in range(len(v_Path))]	
	# fig = plt.figure()
	# ax = fig.add_subplot(111, projection = '3d')
	# ax.scatter(Px_Path, Py_Path, Pz_Path, c = 'm', marker = 'o')
	# plta = plt.gca(projection='3d')
	# plta.plot(Px_Path, Py_Path, Pz_Path)
	plt.plot(x_axis,theta1)
	plt.plot(x_axis,theta2)
	plt.plot(x_axis,theta3)
	plt.plot(x_axis,theta4)
	plt.plot(x_axis,theta5)
	plt.plot(x_axis,theta6)
	# plt.plot(x1_axis,v_Path)

	plt.show()

	return np.array([[theta1],
					[theta2],
					[theta3],
					[theta4],
					[theta5],
					[theta6],
					[Px_Path],
					[Py_Path],
					[Pz_Path],
					[v_Path]])

def MZ07_Trajectory_Planning_StraightLine_Polynomial_Twist_AngleOutput(θ1_origin, θ2_origin,θ3_origin,θ4_origin,θ5_origin,θ6_origin,θ1_target, θ2_target,θ3_target,θ4_target,θ5_target,θ6_target,T=3,t=0.01,steps=10,error=0.001):
	a0 = 0 
	a1 = 0
	a2 = 3/(T*T)
	a3 = -2/(T*T*T)

	M = np.array([[1, 0, 0, 50+330+45],
				  [0,-1, 0,         0],
			 	  [0, 0,-1,345-340-73],
				  [0, 0, 0,         1]])
	
	errorlist = []
	theta1 = []
	theta2 = []
	theta3 = []
	theta4 = []
	theta5 = []
	theta6 = []

	I = np.array([[1,0,0],
				  [0,1,0],
				  [0,0,1]])

	S1 = np.array([[0],[0],[1],[0],[0],[0]])
	S2 = np.array([[0],[-1],[0],[345],[0],[-50]])
	S3 = np.array([[0],[-1],[0],[345],[0],[-330-50]])
	S4 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])
	S5 = np.array([[0],[-1],[0],[345-340],[0],[-50-330-45]])
	S6 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])

	B1 = np.array([[0],[0],[-1],[0],[-50-330-45],[0]])
	B2 = np.array([[0],[1],[0],[340+73],[0],[-330-45]])
	B3 = np.array([[0],[1],[0],[340+73],[0],[-45]])
	B4 = np.array([[0],[0],[1],[0],[0],[0]])
	B5 = np.array([[0],[1],[0],[73],[0],[0]])	
	B6 = np.array([[0],[0],[1],[0],[0],[0]])

	'''设定终点'''
	θ_target = np.array([[θ1_target*np.pi/180],
					  	 [θ2_target*np.pi/180],
						 [θ3_target*np.pi/180],
						 [θ4_target*np.pi/180],
						 [θ5_target*np.pi/180],
						 [θ6_target*np.pi/180]])

	e_Sθ1_target = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ_target[0][0])
	e_Sθ2_target = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ_target[1][0])
	e_Sθ3_target = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ_target[2][0])
	e_Sθ4_target = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ_target[3][0])
	e_Sθ5_target = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],θ_target[4][0])
	e_Sθ6_target = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],θ_target[5][0])

	T12_target = np.dot(e_Sθ1_target, e_Sθ2_target)
	T13_target = np.dot(T12_target, e_Sθ3_target)
	T14_target = np.dot(T13_target, e_Sθ4_target)
	T15_target = np.dot(T14_target, e_Sθ5_target)
	T16_target = np.dot(T15_target, e_Sθ6_target)
	Ts_target = np.dot(T16_target, M) 

	'''设定起点'''
	θ_origin = np.array([[θ1_origin*np.pi/180],
						 [θ2_origin*np.pi/180],
						 [θ3_origin*np.pi/180],
						 [θ4_origin*np.pi/180],
					 	 [θ5_origin*np.pi/180],
						 [θ6_origin*np.pi/180]])

	e_Sθ1_origin = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ_origin[0][0])
	e_Sθ2_origin = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ_origin[1][0])
	e_Sθ3_origin = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ_origin[2][0])
	e_Sθ4_origin = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ_origin[3][0])
	e_Sθ5_origin = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],θ_origin[4][0])
	e_Sθ6_origin = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],θ_origin[5][0])

	T12_origin = np.dot(e_Sθ1_origin, e_Sθ2_origin)
	T13_origin = np.dot(T12_origin, e_Sθ3_origin)
	T14_origin = np.dot(T13_origin, e_Sθ4_origin)
	T15_origin = np.dot(T14_origin, e_Sθ5_origin)
	T16_origin = np.dot(T15_origin, e_Sθ6_origin)
	Ts_origin = np.dot(T16_origin, M) 

	Px_Path = []
	Py_Path = []
	Pz_Path = []
	vx_Path = []
	vy_Path = []
	vz_Path = []
	v_Path = []
	a0 = 0 
	a1 = 0
	a2 = 3/(T*T)
	a3 = -2/(T*T*T)

	for i in np.arange(0,T+t,t):
		if i == 0:

			theta = θ_origin
			theta1.append(theta[0][0]*180/np.pi)
			theta2.append(theta[1][0]*180/np.pi)
			theta3.append(theta[2][0]*180/np.pi)
			theta4.append(theta[3][0]*180/np.pi)
			theta5.append(theta[4][0]*180/np.pi)
			theta6.append(theta[5][0]*180/np.pi)

			Px_Path.append(Ts_origin[0,3])
			Py_Path.append(Ts_origin[1,3])
			Pz_Path.append(Ts_origin[2,3])

			V_Tst = (6*i/(T*T)-6*i*i/(T*T*T))*(Ts_target-Ts_origin)

			vx_Path.append(V_Tst[0,3])
			vy_Path.append(V_Tst[1,3])
			vz_Path.append(V_Tst[2,3])

			v_Path.append(np.sqrt(vx_Path[-1]*vx_Path[-1] + vy_Path[-1]*vy_Path[-1] + vz_Path[-1]*vz_Path[-1]))
			
			i = i + t

		s = a0 + a1*i+a2*i*i+a3*i*i*i
		Tst = Ts_origin + s*(Ts_target-Ts_origin)


		'''逆运动学数值解'''
		for j in range(steps): 

			e_Sθ1 = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],theta[0][0])
			e_Sθ2 = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],theta[1][0])
			e_Sθ3 = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],theta[2][0])
			e_Sθ4 = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],theta[3][0])
			e_Sθ5 = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],theta[4][0])
			e_Sθ6 = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],theta[5][0])

			Ts_b12 = np.dot(e_Sθ1, e_Sθ2)
			Ts_b13 = np.dot(Ts_b12, e_Sθ3)
			Ts_b14 = np.dot(Ts_b13, e_Sθ4)	
			Ts_b15 = np.dot(Ts_b14, e_Sθ5)
			Ts_b16 = np.dot(Ts_b15, e_Sθ6)
			Ts_b = np.dot(Ts_b16,M) 
			#print('Ts_b =', Ts_b)
			Ts_b_Inverse = np.linalg.inv(Ts_b)		
			Tb_s = Ts_b_Inverse 
			Tb_target = np.dot(Tb_s,Tst)
			
			#print('Tb_target =', Tb_target)

			nx = Tb_target[0,0]
			ny = Tb_target[1,0]
			nz = Tb_target[2,0]
			ox = Tb_target[0,1]
			oy = Tb_target[1,1]
			oz = Tb_target[2,1]
			ax = Tb_target[0,2]
			ay = Tb_target[1,2]
			az = Tb_target[2,2]
			px = Tb_target[0,3]
			py = Tb_target[1,3]
			pz = Tb_target[2,3]

			θi = np. arccos((nx + oy + az  - 1)/2)


			ωi1 = (oz-ay)/(2*np.sin(θi))
			ωi2 = (ax-nz)/(2*np.sin(θi))
			ωi3 = (ny-ox)/(2*np.sin(θi))

			ωi = np.random.random(9)
			ωi = ωi.reshape(3, 3)

			ωi[0,0] = 0               
			ωi[0,1] = -ωi3                   
			ωi[0,2] = ωi2        
			ωi[1,0] = ωi3              
			ωi[1,1] = 0                  
			ωi[1,2] = -ωi1    
			ωi[2,0] = -ωi2               
			ωi[2,1] = ωi1                   
			ωi[2,2] = 0  

			#print('ωi =',ωi)
			'''
			Ri_Inverse = I/θi-ωi/2+ np.dot((1/θi-1/(np.tan(θi/2)*2)),np.dot(ωi,ωi))

			'''
			Ri = np.dot(I,θi)+np.dot((1-np.cos(θi)),ωi)+np.dot((θi-np.sin(θi)),np.dot(ωi,ωi))			

			Ri_Inverse = np.linalg.inv(Ri)

			Pi = [[px],
				  [py],
				  [pz]]
			vi = np.dot(Ri_Inverse,Pi)

			vi1 = vi[0][0]
			vi2 = vi[1][0]
			vi3 = vi[2][0]

			Vbi = np.array([[ωi1*θi],
							[ωi2*θi],
							[ωi3*θi],
							[vi1*θi],
							[vi2*θi],
							[vi3*θi]])
		
			AdT_b1 = PoE.AdT(B1[0][0],B1[1][0],B1[2][0],B1[3][0],B1[4][0],B1[5][0],-theta[0][0])
			AdT_b2 = PoE.AdT(B2[0][0],B2[1][0],B2[2][0],B2[3][0],B2[4][0],B2[5][0],-theta[1][0])
			AdT_b3 = PoE.AdT(B3[0][0],B3[1][0],B3[2][0],B3[3][0],B3[4][0],B3[5][0],-theta[2][0])
			AdT_b4 = PoE.AdT(B4[0][0],B4[1][0],B4[2][0],B4[3][0],B4[4][0],B4[5][0],-theta[3][0])
			AdT_b5 = PoE.AdT(B5[0][0],B5[1][0],B5[2][0],B5[3][0],B5[4][0],B5[5][0],-theta[4][0])
			AdT_b6 = PoE.AdT(B6[0][0],B6[1][0],B6[2][0],B6[3][0],B6[4][0],B6[5][0],-theta[5][0])


			AdT_b65 = np.dot(AdT_b6,AdT_b5)		
			AdT_b654= np.dot(AdT_b65,AdT_b4)
			AdT_b6543 = np.dot(AdT_b654,AdT_b3)
			AdT_b65432 = np.dot(AdT_b6543,AdT_b2)


			Jb1 = np.dot(AdT_b65432,B1)
			Jb2 = np.dot(AdT_b6543,B2)
			Jb3 = np.dot(AdT_b654,B3)
			Jb4 = np.dot(AdT_b65,B4)
			Jb5 = np.dot(AdT_b6,B5)
			Jb6 = B6

			Jb = np.random.random(36)
			Jb = Jb.reshape(6, 6)

			Jb[0,0] = Jb1[0][0]
			Jb[1,0] = Jb1[1][0]
			Jb[2,0] = Jb1[2][0]
			Jb[3,0] = Jb1[3][0]
			Jb[4,0] = Jb1[4][0]
			Jb[5,0] = Jb1[5][0]

			Jb[0,1] = Jb2[0][0]
			Jb[1,1] = Jb2[1][0]
			Jb[2,1] = Jb2[2][0]
			Jb[3,1] = Jb2[3][0]
			Jb[4,1] = Jb2[4][0]
			Jb[5,1] = Jb2[5][0]

			Jb[0,2] = Jb3[0][0]
			Jb[1,2] = Jb3[1][0]
			Jb[2,2] = Jb3[2][0]
			Jb[3,2] = Jb3[3][0]
			Jb[4,2] = Jb3[4][0]
			Jb[5,2] = Jb3[5][0]

			Jb[0,3] = Jb4[0][0]
			Jb[1,3] = Jb4[1][0]
			Jb[2,3] = Jb4[2][0]
			Jb[3,3] = Jb4[3][0]
			Jb[4,3] = Jb4[4][0]
			Jb[5,3] = Jb4[5][0]

			Jb[0,4] = Jb5[0][0]
			Jb[1,4] = Jb5[1][0]
			Jb[2,4] = Jb5[2][0]
			Jb[3,4] = Jb5[3][0]
			Jb[4,4] = Jb5[4][0]
			Jb[5,4] = Jb5[5][0]

			Jb[0,5] = Jb6[0][0]
			Jb[1,5] = Jb6[1][0]
			Jb[2,5] = Jb6[2][0]
			Jb[3,5] = Jb6[3][0]
			Jb[4,5] = Jb6[4][0]
			Jb[5,5] = Jb6[5][0]

			'''求矩阵的逆'''
			Jb_Inverse = np.linalg.inv(Jb)

			thetaX = theta
			theta = theta + np.dot(Jb_Inverse,Vbi)
			deltatheta = theta-thetaX

			#errorlist.append(np.linalg.norm(deltatheta))
			e = np.sqrt(np.dot(deltatheta[0][0],deltatheta[0][0]) + np.dot(deltatheta[1][0],deltatheta[1][0]) + np.dot(deltatheta[2][0],deltatheta[2][0]) + np.dot(deltatheta[3][0],deltatheta[3][0])+ np.dot(deltatheta[4][0],deltatheta[4][0])+ np.dot(deltatheta[5][0],deltatheta[5][0]))
			errorlist.append(e)
				
			if errorlist[-1] < error:
				# print('iteration: ', theta*180/np.pi%360)
				print('迭代次数:',j+1)
				theta1.append(theta[0][0]*180/np.pi)
				theta2.append(theta[1][0]*180/np.pi)
				theta3.append(theta[2][0]*180/np.pi)
				theta4.append(theta[3][0]*180/np.pi)
				theta5.append(theta[4][0]*180/np.pi)
				theta6.append(theta[5][0]*180/np.pi)
				break


		Px_Path.append(Tst[0,3])
		Py_Path.append(Tst[1,3])
		Pz_Path.append(Tst[2,3])


		V_Tst = (6*i/(T*T)-6*i*i/(T*T*T))*(Ts_target-Ts_origin)

		vx_Path.append(V_Tst[0,3])
		vy_Path.append(V_Tst[1,3])
		vz_Path.append(V_Tst[2,3])

		v_Path.append(np.sqrt(vx_Path[-1]*vx_Path[-1] + vy_Path[-1]*vy_Path[-1] + vz_Path[-1]*vz_Path[-1]))


	
	x_axis = [i for i in range(len(theta1))]
	x1_axis = [i for i in range(len(v_Path))]	
	plt.figure()
	# plta = plt.gca(projection='3d')
	# plta.plot(Px_Path, Py_Path, Pz_Path)
	plt.plot(x_axis,theta1)
	plt.plot(x_axis,theta2)
	plt.plot(x_axis,theta3)
	plt.plot(x_axis,theta4)
	plt.plot(x_axis,theta5)
	plt.plot(x_axis,theta6)
	# plt.plot(x1_axis,v_Path)

	plt.show()

	return np.array([[theta1],
					[theta2],
					[theta3],
					[theta4],
					[theta5],
					[theta6],
					[Px_Path],
					[Py_Path],
					[Pz_Path],
					[v_Path]])

def MZ07_Trajectory_Planning_StraightLine_TriplePolynomial_Twist(θ1_origin, θ2_origin,θ3_origin,θ4_origin,θ5_origin,θ6_origin,θ1_target, θ2_target,θ3_target,θ4_target,θ5_target,θ6_target,T=3,t=0.01):


	M = np.array([[1, 0, 0, 50+330+45],

				  [0,-1, 0,         0],
			 	  [0, 0,-1,345-340-73],
				  [0, 0, 0,         1]])

	S1 = np.array([[0],[0],[1],[0],[0],[0]])
	S2 = np.array([[0],[-1],[0],[345],[0],[-50]])
	S3 = np.array([[0],[-1],[0],[345],[0],[-330-50]])
	S4 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])
	S5 = np.array([[0],[-1],[0],[345-340],[0],[-50-330-45]])
	S6 = np.array([[0],[0],[-1],[0],[50+330+45],[0]])
	
	'''设定终点'''
	e_Sθ1_target = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ1_target*np.pi/180)
	e_Sθ2_target = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ2_target*np.pi/180)
	e_Sθ3_target = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ3_target*np.pi/180)
	e_Sθ4_target = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ4_target*np.pi/180)
	e_Sθ5_target = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],θ5_target*np.pi/180)
	e_Sθ6_target = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],θ6_target*np.pi/180)

	T12_target = np.dot(e_Sθ1_target, e_Sθ2_target)
	T13_target = np.dot(T12_target, e_Sθ3_target)
	T14_target = np.dot(T13_target, e_Sθ4_target)
	T15_target = np.dot(T14_target, e_Sθ5_target)
	T16_target = np.dot(T15_target, e_Sθ6_target)
	Ts_target = np.dot(T16_target, M) 

	'''设定起点'''
	e_Sθ1_origin = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ1_origin*np.pi/180)
	e_Sθ2_origin = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ2_origin*np.pi/180)
	e_Sθ3_origin = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ3_origin*np.pi/180)
	e_Sθ4_origin = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ4_origin*np.pi/180)
	e_Sθ5_origin = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],θ5_origin*np.pi/180)
	e_Sθ6_origin = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],θ6_origin*np.pi/180)

	T12_origin = np.dot(e_Sθ1_origin, e_Sθ2_origin)
	T13_origin = np.dot(T12_origin, e_Sθ3_origin)
	T14_origin = np.dot(T13_origin, e_Sθ4_origin)
	T15_origin = np.dot(T14_origin, e_Sθ5_origin)
	T16_origin = np.dot(T15_origin, e_Sθ6_origin)
	Ts_origin = np.dot(T16_origin, M) 

	P_target = np.random.random(3)
	P_target = P_target.reshape(3, 1)

	P_target[0][0] = Ts_target[0][3]
	P_target[1][0] = Ts_target[1][3]
	P_target[2][0] = Ts_target[2][3]

	P_origin = np.random.random(3)
	P_origin = P_origin.reshape(3, 1)

	P_origin[0][0] = Ts_origin[0][3]
	P_origin[1][0] = Ts_origin[1][3]
	P_origin[2][0] = Ts_origin[2][3]

	R_target = np.random.random(9)
	R_target = R_target.reshape(3, 3)

	R_target[0,0] = Ts_target[0,0]
	R_target[1,0] = Ts_target[1,0]
	R_target[2,0] = Ts_target[2,0]
	R_target[0,1] = Ts_target[0,1]
	R_target[1,1] = Ts_target[1,1]
	R_target[2,1] = Ts_target[2,1]
	R_target[0,2] = Ts_target[0,2]
	R_target[1,2] = Ts_target[1,2]
	R_target[2,2] = Ts_target[2,2]

	R_origin = np.random.random(9)
	R_origin = R_origin.reshape(3, 3)

	R_origin[0,0] = Ts_origin[0,0]
	R_origin[1,0] = Ts_origin[1,0]
	R_origin[2,0] = Ts_origin[2,0]
	R_origin[0,1] = Ts_origin[0,1]
	R_origin[1,1] = Ts_origin[1,1]
	R_origin[2,1] = Ts_origin[2,1]
	R_origin[0,2] = Ts_origin[0,2]
	R_origin[1,2] = Ts_origin[1,2]
	R_origin[2,2] = Ts_origin[2,2]

	R_origin_Transpose = R_origin.transpose()
	R_ot = np.dot(R_origin_Transpose,R_target)

	nx_ot = R_ot[0,0]
	ny_ot = R_ot[1,0]
	nz_ot = R_ot[2,0]
	ox_ot = R_ot[0,1]
	oy_ot = R_ot[1,1]
	oz_ot = R_ot[2,1]
	ax_ot = R_ot[0,2]
	ay_ot = R_ot[1,2]
	az_ot = R_ot[2,2]

	θ_ot = np.arccos((nx_ot + oy_ot + az_ot  - 1)/2)

	ω_ot1 = (oz_ot-ay_ot)/(2*np.sin(θ_ot))
	ω_ot2 = (ax_ot-nz_ot)/(2*np.sin(θ_ot))
	ω_ot3 = (ny_ot-ox_ot)/(2*np.sin(θ_ot))

	px = []
	py = []
	pz = []
	vx = []
	vy = []
	vz = []
	v = []
	ωθ_st1 = []
	ωθ_st2 = []
	ωθ_st3 = []
	a0 = 0 
	a1 = 0
	a2 = 3/(T*T)
	a3 = -2/(T*T*T)



	for i in np.arange(0,T+t,t):
		s = a0 + a1*i+a2*i*i+a3*i*i*i
		s_Derivative = a1+2*a2*i+3*a3*i*i

		#平移位置
		Pst = P_origin + s*(P_target-P_origin)

		px.append(Pst[0,0])
		py.append(Pst[1,0])
		pz.append(Pst[2,0])

		#平移速度
		V_Pst = s_Derivative*(P_target-P_origin)

		vx.append(V_Pst[0,0])
		vy.append(V_Pst[1,0])
		vz.append(V_Pst[2,0])
		v.append(np.sqrt(vx[-1]*vx[-1] + vy[-1]*vy[-1] + vz[-1]*vz[-1]))

		#旋转位置
		R_ot = PoE.e_ωθ(ω_ot1,ω_ot2,ω_ot3,θ_ot*s)
		Rst = np.dot(R_origin,R_ot)
		
		nx_st = Rst[0,0]
		ny_st = Rst[1,0]
		nz_st = Rst[2,0]
		ox_st = Rst[0,1]
		oy_st = Rst[1,1]
		oz_st = Rst[2,1]
		ax_st = Rst[0,2]
		ay_st = Rst[1,2]
		az_st = Rst[2,2]

		θ_st = np.arccos((nx_st + oy_st + az_st  - 1)/2)

		ω_st1 = (oz_st-ay_st)/(2*np.sin(θ_st))
		ω_st2 = (ax_st-nz_st)/(2*np.sin(θ_st))
		ω_st3 = (ny_st-ox_st)/(2*np.sin(θ_st))

		ωθ_st1.append(ω_st1*θ_st)
		ωθ_st2.append(ω_st2*θ_st)
		ωθ_st3.append(ω_st3*θ_st)

		#旋转速度
		V_R_ot = PoE.e_ωθ(ω_ot1,ω_ot2,ω_ot3,θ_ot*s_Derivative)
		V_Rst = np.dot(R_origin,V_R_ot)
		
		nx_Vst = V_Rst[0,0]
		ny_Vst = V_Rst[1,0]
		nz_Vst = V_Rst[2,0]
		ox_Vst = V_Rst[0,1]
		oy_Vst = V_Rst[1,1]
		oz_Vst = V_Rst[2,1]
		ax_Vst = V_Rst[0,2]
		ay_Vst = V_Rst[1,2]
		az_Vst = V_Rst[2,2]

		θ_Vst = np.arccos((nx_Vst + oy_Vst + az_Vst - 1)/2)

		ω_Vst1 = (oz_Vst-ay_Vst)/(2*np.sin(θ_Vst))
		ω_Vst2 = (ax_Vst-nz_Vst)/(2*np.sin(θ_Vst))
		ω_Vst3 = (ny_Vst-ox_Vst)/(2*np.sin(θ_Vst))

		ωθ_st1.append(ω_Vst1*θ_Vst)
		ωθ_st2.append(ω_Vst2*θ_Vst)
		ωθ_st3.append(ω_Vst3*θ_Vst)

	print(v)
	x_axis = [i for i in range(len(px))]
	plt.figure()
	#plta = plt.gca(projection='3d')
	#plta.plot(px, py, pz)
	plt.plot(x_axis,ωθ_st1)
	plt.plot(x_axis,ωθ_st2)
	plt.plot(x_axis,ωθ_st3)
	plt.show()


EC06_Trajectory_Planning_StraightLine_Polynomial_AngleOutput(10,10,1,10,10,11,10,45,T=3,t=0.01,steps=10000,error=0.001)
