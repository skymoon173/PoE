
# -*- coding: utf-8 -*-
# coding=utf-8
# coding: utf-8

import sympy
import numpy as np
import matplotlib.pyplot as plt

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
	
	def AdT_Matrix(T):

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


	def adv(ω1,ω2,ω3,v1,v2,v3):
		ω = np.array([[0,-ω3,ω2],
					 [ω3,0,-ω1],
					 [-ω2,ω1,0]])

		v = np.array([[0,-v3,v2],
					 [v3,0,-v1],
					 [-v2,v1,0]])

		adv = np.random.random(36)
		adv = adv.reshape(6, 6)

		adv[0,0] = ω[0,0]
		adv[1,0] = ω[1,0]
		adv[2,0] = ω[2,0]
		adv[0,1] = ω[0,1]
		adv[1,1] = ω[1,1]
		adv[2,1] = ω[2,1]
		adv[0,2] = ω[0,2]
		adv[1,2] = ω[1,2]
		adv[2,2] = ω[2,2]

		adv[3,0] = v[0,0]
		adv[4,0] = v[1,0]
		adv[5,0] = v[2,0]
		adv[3,1] = v[0,1]
		adv[4,1] = v[1,1]
		adv[5,1] = v[2,1]
		adv[3,2] = v[0,2]
		adv[4,2] = v[1,2]
		adv[5,2] = v[2,2]

		adv[0,3] = 0
		adv[1,3] = 0
		adv[2,3] = 0
		adv[0,4] = 0
		adv[1,4] = 0
		adv[2,4] = 0
		adv[0,5] = 0
		adv[1,5] = 0
		adv[2,5] = 0

		adv[3,3] = ω[0,0]
		adv[4,3] = ω[1,0]
		adv[5,3] = ω[2,0]
		adv[3,4] = ω[0,1]
		adv[4,4] = ω[1,1]
		adv[5,4] = ω[2,1]
		adv[3,5] = ω[0,2]
		adv[4,5] = ω[1,2]
		adv[5,5] = ω[2,2]		

		return adv

def MZ07_Inverse_Dynamics_PoEModel(θ1,θ2,θ3,θ4,θ5,θ6,θ1_dot,θ2_dot,θ3_dot,θ4_dot,θ5_dot,θ6_dot,θ1_doubleDot,θ2_doubleDot,θ3_doubleDot,θ4_doubleDot,θ5_doubleDot,θ6_doubleDot):

	S1 = np.array([[0],[0],[1],[0],[0],[0]])
	S2 = np.array([[0],[-1],[0],[0.345],[0],[-0.05]])
	S3 = np.array([[0],[-1],[0],[0.345],[0],[-0.33-0.05]])
	S4 = np.array([[0],[0],[-1],[0],[0.05+0.33+0.045],[0]])
	S5 = np.array([[0],[-1],[0],[0.345-0.34],[0],[-0.05-0.33-0.045]])
	S6 = np.array([[0],[0],[-1],[0],[0.05+0.33+0.045],[0]])

	#定义各连杆质心的位置，大小以及惯性积

	M01 = np.array([[1,0,0, 0.01],
					[0,1,0,    0],
					[0,0,1, 0.28],
					[0,0,0,    1]])

	M02 = np.array([[1,0,0, 0.21],
					[0,1,0,    0],
					[0,0,1, 0.36],
					[0,0,0,    1]])	

	M03 = np.array([[1,0,0, 0.4],
					[0,1,0,   0],
					[0,0,1,0.35],
					[0,0,0,   1]])

	M04 = np.array([[1,0,0,0.42],
					[0,1,0,   0],
					[0,0,1,0.15],
					[0,0,0,   1]])	

	M05 = np.array([[1,0,0, 0.42],
					[0,1,0,    0],
					[0,0,1,-0.02],
					[0,0,0,    1]])

	M06 = np.array([[1,0,0, 0.43],
					[0,1,0,    0],
					[0,0,1,-0.07],
					[0,0,0,    1]])

	m6 = 0.23635
	Ixx6 = 0
	Ixy6 = 0
	Ixz6 = 0
	Iyy6 = 0
	Iyz6 = 0
	Izz6 = 0.00016748

	m5 = 1.0645
	Ixx5 = 0.0014536
	Ixy5 = 0
	Ixz5 = 0
	Iyy5 = 0.0013263
	Iyz5 = 0.0001127
	Izz5 = 0.0011216

	m4 = 6.8836
	Ixx4 = 0.040912
	Ixy4 = 0
	Ixz4 = 0
	Iyy4 = 0.036854
	Iyz4 = 0.00049331
	Izz4 = 0.015137

	m3 = 7.8998
	Ixx3 = 0.02937
	Ixy3 = 0
	Ixz3 = 0.0024455
	Iyy3 = 0.037441
	Iyz3 = 0
	Izz3 = 0.021547

	m2 = 15.589
	Ixx2 = 0.074369
	Ixy2 = 0
	Ixz2 = 0
	Iyy2 = 0.2114
	Iyz2 = 0
	Izz2 = 0.25542

	m1 = 11.913
	Ixx1 = 0.096913
	Ixy1 = 0
	Ixz1 = -0.01294
	Iyy1 = 0.063825
	Iyz1 = 0
	Izz1 = 0.081833

	M01_Inverse = np.linalg.inv(M01)
	M02_Inverse = np.linalg.inv(M02)
	M03_Inverse = np.linalg.inv(M03)
	M04_Inverse = np.linalg.inv(M04)
	M05_Inverse = np.linalg.inv(M05)
	M06_Inverse = np.linalg.inv(M06) 

	M10 = M01_Inverse
	M21 = np.dot(M02_Inverse,M01)
	M32 = np.dot(M03_Inverse,M02)
	M43 = np.dot(M04_Inverse,M03)
	M54 = np.dot(M05_Inverse,M04)
	M65 = np.dot(M06_Inverse,M05)

	AdT_M01_Inverse = PoE.AdT_Matrix(M01_Inverse)
	AdT_M02_Inverse = PoE.AdT_Matrix(M02_Inverse)
	AdT_M03_Inverse = PoE.AdT_Matrix(M03_Inverse)
	AdT_M04_Inverse = PoE.AdT_Matrix(M04_Inverse)
	AdT_M05_Inverse = PoE.AdT_Matrix(M05_Inverse)
	AdT_M06_Inverse = PoE.AdT_Matrix(M06_Inverse)

	#计算质心坐标系下关节的旋量
	A1 = np.dot(AdT_M01_Inverse,S1)
	A2 = np.dot(AdT_M02_Inverse,S2)
	A3 = np.dot(AdT_M03_Inverse,S3)
	A4 = np.dot(AdT_M04_Inverse,S4)
	A5 = np.dot(AdT_M05_Inverse,S5)
	A6 = np.dot(AdT_M06_Inverse,S6)

	#计算后一个连杆质心坐标系在前一个连杆质心坐标系下的位形（串联式连接）
	e_Aθ1_minus = PoE.e_Sθ(A1[0][0],A1[1][0],A1[2][0],A1[3][0],A1[4][0],A1[5][0],-θ1*np.pi/180)
	e_Aθ2_minus = PoE.e_Sθ(A2[0][0],A2[1][0],A2[2][0],A2[3][0],A2[4][0],A2[5][0],-θ2*np.pi/180)
	e_Aθ3_minus = PoE.e_Sθ(A3[0][0],A3[1][0],A3[2][0],A3[3][0],A3[4][0],A3[5][0],-θ3*np.pi/180)
	e_Aθ4_minus = PoE.e_Sθ(A4[0][0],A4[1][0],A4[2][0],A4[3][0],A4[4][0],A4[5][0],-θ4*np.pi/180)
	e_Aθ5_minus = PoE.e_Sθ(A5[0][0],A5[1][0],A5[2][0],A5[3][0],A5[4][0],A5[5][0],-θ5*np.pi/180)
	e_Aθ6_minus = PoE.e_Sθ(A6[0][0],A6[1][0],A6[2][0],A6[3][0],A6[4][0],A6[5][0],-θ6*np.pi/180)

	T10 = np.dot(e_Aθ1_minus,M10)
	T21 = np.dot(e_Aθ2_minus,M21)
	T32 = np.dot(e_Aθ3_minus,M32)
	T43 = np.dot(e_Aθ4_minus,M43)
	T54 = np.dot(e_Aθ5_minus,M54)
	T65 = np.dot(e_Aθ6_minus,M65)
	
	#计算各连杆运动旋量
	Ad_T10 = PoE.AdT_Matrix(T10)
	Ad_T21 = PoE.AdT_Matrix(T21)
	Ad_T32 = PoE.AdT_Matrix(T32)
	Ad_T43 = PoE.AdT_Matrix(T43)
	Ad_T54 = PoE.AdT_Matrix(T54)
	Ad_T65 = PoE.AdT_Matrix(T65)

	V0 = np.array([[0],[0],[0],[0],[0],[0]])
	V1 = np.dot(Ad_T10,V0) + A1*θ1_dot
	V2 = np.dot(Ad_T21,V1) + A2*θ2_dot
	V3 = np.dot(Ad_T32,V2) + A3*θ3_dot
	V4 = np.dot(Ad_T43,V3) + A4*θ4_dot
	V5 = np.dot(Ad_T54,V4) + A5*θ5_dot
	V6 = np.dot(Ad_T65,V5) + A6*θ6_dot

	#计算各连杆运动旋量关于时间的微分
	adv1 = PoE.adv(A1[0][0],A1[1][0],A1[2][0],A1[3][0],A1[4][0],A1[5][0])
	adv2 = PoE.adv(A2[0][0],A2[1][0],A2[2][0],A2[3][0],A2[4][0],A2[5][0])
	adv3 = PoE.adv(A3[0][0],A3[1][0],A3[2][0],A3[3][0],A3[4][0],A3[5][0])
	adv4 = PoE.adv(A4[0][0],A4[1][0],A4[2][0],A4[3][0],A4[4][0],A4[5][0])
	adv5 = PoE.adv(A5[0][0],A5[1][0],A5[2][0],A5[3][0],A5[4][0],A5[5][0])
	adv6 = PoE.adv(A6[0][0],A6[1][0],A6[2][0],A6[3][0],A6[4][0],A6[5][0])

	g = 9.8
	V0_dot = np.array([[0],[0],[0],[0],[0],[-g]])
	V1_dot = np.dot(Ad_T10,V0_dot) + np.dot(adv1,A1)*θ1_dot + A1*θ1_doubleDot
	V2_dot = np.dot(Ad_T21,V1_dot) + np.dot(adv2,A2)*θ2_dot + A2*θ2_doubleDot
	V3_dot = np.dot(Ad_T32,V2_dot) + np.dot(adv3,A3)*θ3_dot + A3*θ3_doubleDot
	V4_dot = np.dot(Ad_T43,V3_dot) + np.dot(adv4,A4)*θ4_dot + A4*θ4_doubleDot
	V5_dot = np.dot(Ad_T54,V4_dot) + np.dot(adv5,A5)*θ5_dot + A5*θ5_doubleDot
	V6_dot = np.dot(Ad_T65,V5_dot) + np.dot(adv6,A6)*θ6_dot + A6*θ6_doubleDot
	
	#计算产生各连杆运动旋量的力旋量
	Ad_T65_Transpose = np.transpose(Ad_T65)
	Ad_T54_Transpose = np.transpose(Ad_T54)
	Ad_T43_Transpose = np.transpose(Ad_T43)
	Ad_T32_Transpose = np.transpose(Ad_T32)
	Ad_T21_Transpose = np.transpose(Ad_T21)

	adv1_Transpose = np.transpose(adv1)
	adv2_Transpose = np.transpose(adv2)
	adv3_Transpose = np.transpose(adv3)
	adv4_Transpose = np.transpose(adv4)
	adv5_Transpose = np.transpose(adv5)
	adv6_Transpose = np.transpose(adv6)

	G6 = np.array([[Ixx6, Ixy6, Ixz6,0,0,0],
				   [Ixy6, Iyy6, Iyz6,0,0,0],
				   [Ixz6, Iyz6, Izz6,0,0,0],
				   [0,0,0,m6,m6,m6],
				   [0,0,0,m6,m6,m6],
				   [0,0,0,m6,m6,m6]]) 

	G5 = np.array([[Ixx5, Ixy5, Ixz5,0,0,0],
				   [Ixy5, Iyy5, Iyz5,0,0,0],
				   [Ixz5, Iyz5, Izz5,0,0,0],
				   [0,0,0,m5,m5,m5],
				   [0,0,0,m5,m5,m5],
				   [0,0,0,m5,m5,m5]]) 	

	G4 = np.array([[Ixx4, Ixy4, Ixz4,0,0,0],
				   [Ixy4, Iyy4, Iyz4,0,0,0],
				   [Ixz4, Iyz4, Izz4,0,0,0],
				   [0,0,0,m4,m4,m4],
				   [0,0,0,m4,m4,m4],
				   [0,0,0,m4,m4,m4]]) 

	G3 = np.array([[Ixx3, Ixy3, Ixz3,0,0,0],
				   [Ixy3, Iyy3, Iyz3,0,0,0],
				   [Ixz3, Iyz3, Izz3,0,0,0],
				   [0,0,0,m3,m3,m3],
				   [0,0,0,m3,m3,m3],
				   [0,0,0,m3,m3,m3]]) 

	G2 = np.array([[Ixx2, Ixy2, Ixz2,0,0,0],
				   [Ixy2, Iyy2, Iyz2,0,0,0],
				   [Ixz2, Iyz2, Izz2,0,0,0],
				   [0,0,0,m2,m2,m2],
				   [0,0,0,m2,m2,m2],
				   [0,0,0,m2,m2,m2]]) 

	G1 = np.array([[Ixx1, Ixy1, Ixz1,0,0,0],
				   [Ixy1, Iyy1, Iyz1,0,0,0],
				   [Ixz1, Iyz1, Izz1,0,0,0],
				   [0,0,0,m1,m1,m1],
				   [0,0,0,m1,m1,m1],
				   [0,0,0,m1,m1,m1]]) 
				   

	F6 = np.dot(G6,V6_dot)-np.dot(adv6_Transpose,np.dot(G6,V6))
	F5 = np.dot(Ad_T65_Transpose,F6) + np.dot(G5,V5_dot)-np.dot(adv5_Transpose,np.dot(G5,V5))
	F4 = np.dot(Ad_T54_Transpose,F5) + np.dot(G4,V4_dot)-np.dot(adv4_Transpose,np.dot(G4,V4))
	F3 = np.dot(Ad_T43_Transpose,F4) + np.dot(G3,V3_dot)-np.dot(adv3_Transpose,np.dot(G3,V3))
	F2 = np.dot(Ad_T32_Transpose,F3) + np.dot(G2,V2_dot)-np.dot(adv2_Transpose,np.dot(G2,V2))
	F1 = np.dot(Ad_T21_Transpose,F2) + np.dot(G1,V1_dot)-np.dot(adv1_Transpose,np.dot(G1,V1))

	#计算各轴力矩
	F6_Transpose = np.transpose(F6)
	F5_Transpose = np.transpose(F5)
	F4_Transpose = np.transpose(F4)
	F3_Transpose = np.transpose(F3)
	F2_Transpose = np.transpose(F2)
	F1_Transpose = np.transpose(F1)

	τ6 = np.dot(F6_Transpose,A6)
	τ5 = np.dot(F5_Transpose,A5)
	τ4 = np.dot(F4_Transpose,A4)
	τ3 = np.dot(F3_Transpose,A3)
	τ2 = np.dot(F2_Transpose,A2)
	τ1 = np.dot(F1_Transpose,A1)

	τ = np.array([τ1,τ2,τ3,τ4,τ5,τ6])

	τ_Value = np.array([τ1[0,0],τ2[0,0],τ3[0,0],τ4[0,0],τ5[0,0],τ6[0,0]])
	print('τ_Value = ', τ_Value)

	return τ_Value

MZ07_Inverse_Dynamics_PoEModel(0,90,0,0,-90,0,2,2,2,2,2,2,3,3,3,3,3,3)


# def main():
# 	MZ07_Inverse_Dynamics_PoEModel()

# if __name__ == '__main__':
#     main()

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
	# print('T', T)
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


