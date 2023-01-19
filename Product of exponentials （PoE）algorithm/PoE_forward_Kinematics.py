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
	def e_Sθ(ω1,ω2,ω3,v1,v2,v3,θ):
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
	
	S1_T = S1.T[0]
	S2_T = S2.T[0]
	S3_T = S3.T[0]
	S4_T = S4.T[0]
	S5_T = S5.T[0]
	S6_T = S6.T[0]
	'''
	e_Sθ1 = PoE.e_Sθ(0,0,1,0,0,0,θ1*np.pi/180)
	e_Sθ2 = PoE.e_Sθ(0,-1,0,345,0,-50,θ2*np.pi/180)
	e_Sθ3 = PoE.e_Sθ(0,-1,0,345,0,-330-50,θ3*np.pi/180)
	e_Sθ4 = PoE.e_Sθ(0,0,-1,0,50+330+45,0,θ4*np.pi/180)
	e_Sθ5 = PoE.e_Sθ(0,-1,0,345-340,0,-50-330-45,θ5*np.pi/180)
	e_Sθ6 = PoE.e_Sθ(0,0,-1,0,50+330+45,0,θ6*np.pi/180)
	'''
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

	AdT1 = PoE.AdT(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],θ1*np.pi/180)
	AdT2 = PoE.AdT(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],θ2*np.pi/180)
	AdT3 = PoE.AdT(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],θ3*np.pi/180)
	AdT4 = PoE.AdT(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],θ4*np.pi/180)
	AdT5 = PoE.AdT(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],θ5*np.pi/180)
	AdT6 = PoE.AdT(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],θ6*np.pi/180)

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

	'''求矩阵的逆'''
	Js_Inverse = np.linalg.inv(Js)
	'''求矩阵的伪逆(瘦高型)'''	
	Js_T = Js.T	
	Js_pseudoInverse = np.dot(np.linalg.inv(np.dot(Js_T,Js)),Js_T)

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

def EC06_Numerical_Inverse_Kinematics_PoEModel(X,Y,Z,θ,steps=10,theta0_1=5,theta0_2=5,theta0_3=150,theta0_4=5):

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
	
	error = 1

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

	'''根据T04转换矩阵求运动旋量的V'''
	θ = θ*np.pi/180
	v1 = (np.cos(θ)*Y - X*np.sin(θ) - Y)/(np.cos(θ)*np.cos(θ) - 2*np.cos(θ) - np.sin(θ)*np.sin(θ) + 1)
	v2 = (np.cos(θ)*X - X - Y*np.sin(θ))/(np.cos(θ)*np.cos(θ) - 2*np.cos(θ) - np.sin(θ)*np.sin(θ) + 1)
	v3 =  Z/θ

	'''目标位置'''
	xd_target = np.array([[0],[0],[1],[v1],[v2],[v3]])	
	print('xd_target = ',xd_target)
	
	fw1 = open(r"C:\Users\eastw\Desktop\development\Product of exponentials （PoE）algorithm\1.txt", 'w')
	fw2 = open(r"C:\Users\eastw\Desktop\development\Product of exponentials （PoE）algorithm\2.txt", 'w')
	fw3 = open(r"C:\Users\eastw\Desktop\development\Product of exponentials （PoE）algorithm\3.txt", 'w')
	fw4 = open(r"C:\Users\eastw\Desktop\development\Product of exponentials （PoE）algorithm\4.txt", 'w')

	for i in range(steps): 

		e_Sθ1 = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],theta[0][0])
		e_Sθ2 = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],theta[1][0])
		e_Sθ3 = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],theta[2][0])
		e_Sθ4 = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],theta[3][0])

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

		V1 = (oy*py + px*ox - py)/(oy*oy - 2*oy - ox*ox + 1)
		V2 = (ox*py+oy*px-px)/(oy*oy-ox*ox-oy*2+1)
		V3 = pz/np.arctan(-ox/oy)

		Ftheta = [[0],
				  [0],
				  [1],
				  [V1],
				  [V2],
				  [V3]]

		print('Ftheta = ',Ftheta)

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

		J = np.array([[0,0,0,0],
					  [0,0,0,0],
					  [1,1,0,1],
					  [0,325*np.sin(theta[0][0]),0, 325*np.sin(theta[0][0])+275*np.sin(theta[0][0]+theta[1][0])],
					  [0,-325*np.cos(theta[0][0]),0,-325*np.cos(theta[0][0])-275*np.cos(theta[0][0]+theta[1][0])],
					  [0,0,1,0]])

		print('J = ', J)
		print('Js = ', Js)
		'''求矩阵的伪逆(瘦高型)'''	
		Js_T = Js.T	
		Js_pseudoInverse = np.dot(np.linalg.inv(np.dot(Js_T,Js)),Js_T)

		print('设定角度 = ',theta[0][0]*180/np.pi,theta[1][0]*180/np.pi,theta[2][0],theta[3][0]*180/np.pi)
		thetaX = theta%np.pi
		theta = theta + np.dot(Js_pseudoInverse,xd_target-Ftheta)
		theta = theta%np.pi
		deltatheta = theta*180/np.pi-thetaX*180/np.pi

		delta = xd_target-Ftheta

		e1 = np.sqrt(np.dot(delta[0][0],delta[0][0]) + np.dot(delta[1][0] ,delta[1][0] ) + np.dot(delta[2][0] ,delta[2][0] ) + np.dot(delta[3][0] ,delta[3][0] ) + np.dot(delta[4][0] ,delta[4][0] ) + np.dot(delta[5][0] ,delta[5][0] ))
		e2 = np.sqrt(np.dot(deltatheta[0][0],deltatheta[0][0]) + np.dot(deltatheta[1][0] ,deltatheta[1][0] ) + np.dot(deltatheta[2][0] ,deltatheta[2][0] ) + np.dot(deltatheta[3][0] ,deltatheta[3][0]))
		e = e1 +e2
		'''

		# errorlist.append(np.linalg.norm(deltatheta))

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
	plt.plot(x_axis, errorlist)
	plt.plot(x_axis, theta1)	
	plt.plot(x_axis, theta2)
	plt.plot(x_axis, theta3)
	plt.plot(x_axis, theta4)

	plt.show()
	
	return [theta[0][0]*180/np.pi,theta[1][0]*180/np.pi,theta[2][0],theta[3][0]*180/np.pi]

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
	
	r = np.array([[r11,r12,r13],
					[r21,r22,r23],
					[r31,r32,r33]])
	
	θ = np.arccos((r11+r22+r33-1)/2)

	ω1 = (r32-r23)/(2*np.sin(θ))
	ω2 = (r13-r31)/(2*np.sin(θ))
	ω3 = (r21-r12)/(2*np.sin(θ))
	'''
	θ = np.arccos((8.57311628e-01-9.79462371e-01-8.66927689e-01-1)/2)

	ω1 = (-1.44424400e-01-5.66980262e-02)/(2*np.sin(θ))
	ω2 = (4.95198663e-01-4.77051333e-01)/(2*np.sin(θ))
	ω3 = (-1.93491078e-01+1.40694904e-01)/(2*np.sin(θ))
	'''

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

	R = np.dot(I,θ)+np.dot((1-np.cos(θ)),ω)+np.dot((θ-np.sin(θ)),np.dot(ω,ω))
	R_Inverse = np.linalg.inv(R)
	P = [[X],
		 [Y],
		 [Z]]
	V = np.dot(R_Inverse,P)

	v1 = V[0][0]
	v2 = V[1][0]
	v3 = V[2][0]

	'''目标位置'''
	xd_target = np.array([[ω1],[ω2],[ω3],[v1],[v2],[v3]])	

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

	for i in range(steps): 

		e_Sθ1 = PoE.e_Sθ(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],theta[0][0])
		e_Sθ2 = PoE.e_Sθ(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],theta[1][0])
		e_Sθ3 = PoE.e_Sθ(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],theta[2][0])
		e_Sθ4 = PoE.e_Sθ(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],theta[3][0])
		e_Sθ5 = PoE.e_Sθ(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],theta[4][0])
		e_Sθ6 = PoE.e_Sθ(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],theta[5][0])

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

		Ri = np.dot(I,θi)+np.dot((1-np.cos(θi)),ωi)+np.dot((θ-np.sin(θi)),np.dot(ωi,ωi))
		Ri_Inverse = np.linalg.inv(Ri)
		Pi = [[px],
			 [py],
			 [pz]]
		Vi = np.dot(Ri_Inverse,Pi)

		vi1 = Vi[0][0]
		vi2 = Vi[1][0]
		vi3 = Vi[2][0]

		Ftheta = np.array([[ωi1],[ωi2],[ωi3],[vi1],[vi2],[vi3]])
	
		AdT1 = PoE.AdT(S1[0][0],S1[1][0],S1[2][0],S1[3][0],S1[4][0],S1[5][0],theta[0][0])
		AdT2 = PoE.AdT(S2[0][0],S2[1][0],S2[2][0],S2[3][0],S2[4][0],S2[5][0],theta[1][0])
		AdT3 = PoE.AdT(S3[0][0],S3[1][0],S3[2][0],S3[3][0],S3[4][0],S3[5][0],theta[2][0])
		AdT4 = PoE.AdT(S4[0][0],S4[1][0],S4[2][0],S4[3][0],S4[4][0],S4[5][0],theta[3][0])
		AdT5 = PoE.AdT(S5[0][0],S5[1][0],S5[2][0],S5[3][0],S5[4][0],S5[5][0],theta[4][0])
		AdT6 = PoE.AdT(S6[0][0],S6[1][0],S6[2][0],S6[3][0],S6[4][0],S6[5][0],theta[5][0])


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

		'''求矩阵的逆'''
		Js_Inverse = np.linalg.inv(Js)

		'''求矩阵的伪逆(瘦高型)'''	
		Js_T = Js.T	
		Js_pseudoInverse = np.dot(np.linalg.inv(np.dot(Js_T,Js)),Js_T)

		thetaX = theta%np.pi
		theta = theta + np.dot(Js_pseudoInverse,xd_target-Ftheta)
		theta = theta%np.pi
		deltatheta = theta-thetaX
		delta = xd_target-Ftheta

		#errorlist.append(np.linalg.norm(deltatheta))
		e1 = np.sqrt(np.dot(deltatheta[0][0],deltatheta[0][0]) + np.dot(deltatheta[1][0],deltatheta[1][0]) + np.dot(deltatheta[2][0],deltatheta[2][0]) + np.dot(deltatheta[3][0],deltatheta[3][0])+ np.dot(deltatheta[4][0],deltatheta[4][0])+ np.dot(deltatheta[5][0],deltatheta[5][0]))
		e = np.sqrt(np.dot(delta[0][0],delta [0][0]) + np.dot(delta[1][0],delta[1][0]) + np.dot(delta[2][0],delta[2][0]) + np.dot(delta[3][0],delta[3][0])+ np.dot(delta[4][0],delta[4][0])+ np.dot(delta[5][0],delta[5][0]))

		errorlist.append(e)
			
		theta1.append(theta[0][0]*180/np.pi)
		theta2.append(theta[1][0]*180/np.pi)
		theta3.append(theta[2][0]*180/np.pi)
		theta4.append(theta[3][0]*180/np.pi)
		theta5.append(theta[4][0]*180/np.pi)
		theta6.append(theta[5][0]*180/np.pi)

		if errorlist[-1] < error:
			print('iteration: ', theta)
			
			break
	x_axis = [i for i in range(len(errorlist))]
	
	plt.figure()
	plt.plot(x_axis, errorlist)
	#plt.plot(x_axis, theta1)	
	#plt.plot(x_axis, theta2)
	#plt.plot(x_axis, theta3)
	#plt.plot(x_axis, theta4)
	#plt.plot(x_axis, theta5)
	#plt.plot(x_axis, theta6)
	plt.show()

	return [theta[0][0]*180/np.pi,theta[1][0]*180/np.pi,theta[2][0]*180/np.pi,theta[3][0]*180/np.pi,theta[4][0]*180/np.pi,theta[5][0]*180/np.pi]

