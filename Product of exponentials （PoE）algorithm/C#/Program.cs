using System;
using MathNet.Numerics.LinearAlgebra;
using System.Collections;
using System.Collections.Generic;

namespace PoE_Kinematics
{
    class PoE_Algorithm
    {
        static void Main(string[] args)
        {

			PoE_Algorithm n = new PoE_Algorithm();
			Matrix<float> MZ07 = n. MZ07_Trajectory_Planning_StraightLine_TriplePolynomial_AngleOutput(10,20,30,40,50,60,60,50,40,30,20,10);


		}

        public Matrix<double> e_ωθ(double ω1, double ω2, double ω3, double θ)
        {
			var Mat = Matrix<double>.Build;
			double[,] x1 = {{ 1, 0, 0 },
              			    { 0, 1, 0 },
						    { 0, 0, 1 }};
			Matrix<double> I = Mat.DenseOfArray(x1);	
			double[,] x2 = {{ 0, -ω3, ω2 },
              			    { ω3, 0, -ω1 },
						    { -ω2, ω1, 0 }};
			
			Matrix<double> ω = Mat.DenseOfArray(x2);	         
			Matrix<double> R = I + Math.Sin(θ)*ω + (1-Math.Cos(θ))*ω * ω;
			return R;
		}
		public Matrix<double> e_Sθ(double ω1, double ω2, double ω3, double v1, double v2, double  v3, double θ)
		{
			var Mat = Matrix<double>.Build;
			double[,] x1 = {{ 1, 0, 0 },
              			    { 0, 1, 0 },
						    { 0, 0, 1 }};
			Matrix<double> I = Mat.DenseOfArray(x1);	

			double[,] x2 = {{ 0, -ω3, ω2 },
              			    { ω3, 0, -ω1 },
						    { -ω2, ω1, 0 }};			
			Matrix<double> ω = Mat.DenseOfArray(x2);			

			double[,] x3 = {{ v1 },
              			    { v2 },
						    { v3 }};										
			Matrix<double> v = Mat.DenseOfArray(x3);			   

			Matrix<double> e_Sθ = Mat.Dense(4,4);
			PoE_Algorithm PoE = new PoE_Algorithm();
			Matrix<double> R = PoE.e_ωθ(ω1,ω2,ω3,θ);
			e_Sθ[0,0] = R[0,0];
			e_Sθ[1,0] = R[1,0];
			e_Sθ[2,0] = R[2,0];
			e_Sθ[3,0] = 0;

			e_Sθ[0,1] = R[0,1];
			e_Sθ[1,1] = R[1,1];
			e_Sθ[2,1] = R[2,1];
			e_Sθ[3,1] = 0;

			e_Sθ[0,2] = R[0,2];
			e_Sθ[1,2] = R[1,2];
			e_Sθ[2,2] = R[2,2];
			e_Sθ[3,2] = 0;

			Matrix<double> V = (I*θ+(1-Math.Cos(θ))*ω+(θ-Math.Sin(θ))*ω*ω)*v;
			e_Sθ[0,3] = V[0,0];
			e_Sθ[1,3] = V[1,0];
			e_Sθ[2,3] = V[2,0];
			e_Sθ[3,3] = 1;
		
			return e_Sθ;
		}

		public Matrix<double> AdT(double ω1, double ω2, double ω3, double v1, double v2, double  v3, double θ)	
		{
			PoE_Algorithm PoE = new PoE_Algorithm();
			Matrix<double> T = PoE.e_Sθ(ω1,ω2,ω3,v1,v2,v3,θ);
			var Mat = Matrix<double>.Build;
			Matrix<double> R = Mat.Dense(3,3);
			R[0,0] = T[0,0];
			R[1,0] = T[1,0];
			R[2,0] = T[2,0];

			R[0,1] = T[0,1];
			R[1,1] = T[1,1];
			R[2,1] = T[2,1];

			R[0,2] = T[0,2];
			R[1,2] = T[1,2];
			R[2,2] = T[2,2];

			Matrix<double> p = Mat.Dense(3,3);
			p[0,0] = T[0,3];
			p[1,0] = T[1,3];
			p[2,0] = T[2,3];

			Matrix<double> P = Mat.Dense(3,3);
			P[0,0] = 0;   
			P[0,1] = -p[2,0];          
			P[0,2] = p[1,0];        
			P[1,0] = p[2,0];             
			P[1,1] = 0;                  
			P[1,2] = -p[0,0];   
			P[2,0] = -p[1,0];             
			P[2,1] = p[0,0];                  
			P[2,2] = 0;  
			Matrix<double> PR = P*R;

			Matrix<double> AdT = Mat.Dense(6,6);
			AdT[0,0] = R[0,0];
			AdT[1,0] = R[1,0];
			AdT[2,0] = R[2,0];
			AdT[0,1] = R[0,1];
			AdT[1,1] = R[1,1];
			AdT[2,1] = R[2,1];
			AdT[0,2] = R[0,2];
			AdT[1,2] = R[1,2];
			AdT[2,2] = R[2,2];

			AdT[3,0] = PR[0,0];
			AdT[4,0] = PR[1,0];
			AdT[5,0] = PR[2,0];
			AdT[3,1] = PR[0,1];
			AdT[4,1] = PR[1,1];
			AdT[5,1] = PR[2,1];
			AdT[3,2] = PR[0,2];
			AdT[4,2] = PR[1,2];
			AdT[5,2] = PR[2,2];

			AdT[0,3] = 0;
			AdT[1,3] = 0;
			AdT[2,3] = 0;
			AdT[0,4] = 0;
			AdT[1,4] = 0;
			AdT[2,4] = 0;
			AdT[0,5] = 0;
			AdT[1,5] = 0;
			AdT[2,5] = 0;

			AdT[3,3] = R[0,0];
			AdT[4,3] = R[1,0];
			AdT[5,3] = R[2,0];
			AdT[3,4] = R[0,1];
			AdT[4,4] = R[1,1];
			AdT[5,4] = R[2,1];
			AdT[3,5] = R[0,2];
			AdT[4,5] = R[1,2];
			AdT[5,5] = R[2,2];

			return AdT;
		}
		public Vector<double> MZ07_Forward_Kinematics_PoEModel(double θ1, double θ2, double θ3, double θ4, double θ5, double θ6)
		{
			var Mat = Matrix<double>.Build;
			var Vec = Vector<double>.Build;
			Matrix<double> M = Mat.Dense(4,4);					
			M[0,0] = 1;
			M[1,0] = 0;
			M[2,0] = 0;
			M[3,0] = 0;

			M[0,1] = 0;
			M[1,1] = -1;
			M[2,1] = 0;
			M[3,1] = 0;

			M[0,2] = 0;
			M[1,2] = 0;
			M[2,2] = -1;
			M[3,2] = 0;

			M[0,3] = 50+330+45;
			M[1,3] = 0;
			M[2,3] = 345-340-73;
			M[3,3] = 1;

			double[] x1 = { 0 ,
              			    0 ,
						    1 ,
							0 ,
              			    0 ,
							0 };
			Vector<double> S1 = Vec.DenseOfArray(x1);	

			double[] x2 = { 0 ,
              			   -1,
						    0 ,
							345 ,
              			    0 ,
						   -50 };
			Vector<double> S2 = Vec.DenseOfArray(x2);	

			double[] x3 = { 0 ,
              			   -1 ,
							0 ,
							345 ,
              			    0 ,
						   -330-50 };
			Vector<double> S3 = Vec.DenseOfArray(x3);	

			double[] x4 = { 0 ,
              			    0 ,
						   -1 ,
							0 ,
              			    50+330+45 ,
							0 };
			Vector<double> S4 = Vec.DenseOfArray(x4);	

			double[] x5 = { 0 ,
              			   -1 ,
						    0 ,
							345-340 ,
              			    0 ,
						   -50-330-45 };
			Vector<double> S5 = Vec.DenseOfArray(x5);	

			double[] x6 = { 0 ,
              			    0 ,
						   -1 ,
							0 ,
              			    50+330+45 ,
							0 };
			Vector<double> S6 = Vec.DenseOfArray(x6);	

			PoE_Algorithm PoE = new PoE_Algorithm();

			Matrix<double> e_Sθ1 = PoE.e_Sθ(S1[0],S1[1],S1[2],S1[3],S1[4],S1[5],θ1*Math.PI/180);
			Matrix<double> e_Sθ2 = PoE.e_Sθ(S2[0],S2[1],S2[2],S2[3],S2[4],S2[5],θ2*Math.PI/180);
			Matrix<double> e_Sθ3 = PoE.e_Sθ(S3[0],S3[1],S3[2],S3[3],S3[4],S3[5],θ3*Math.PI/180);
			Matrix<double> e_Sθ4 = PoE.e_Sθ(S4[0],S4[1],S4[2],S4[3],S4[4],S4[5],θ4*Math.PI/180);
			Matrix<double> e_Sθ5 = PoE.e_Sθ(S5[0],S5[1],S5[2],S5[3],S5[4],S5[5],θ5*Math.PI/180);
			Matrix<double> e_Sθ6 = PoE.e_Sθ(S6[0],S6[1],S6[2],S6[3],S6[4],S6[5],θ6*Math.PI/180);

			Matrix<double> T12 = e_Sθ1*e_Sθ2;
			Matrix<double> T13 = T12*e_Sθ3;
			Matrix<double> T14 = T13*e_Sθ4;
			Matrix<double> T15 = T14*e_Sθ5;
			Matrix<double> T16 = T15*e_Sθ6;
			Matrix<double> T = T16*M;

			double nx = T[0,0];
			double ny = T[1,0];
			double nz = T[2,0];
			double ox = T[0,1];
			double oy = T[1,1];
			double oz = T[2,1];
			double ax = T[0,2];
			double ay = T[1,2];
			double az = T[2,2];
			double px = T[0,3];
			double py = T[1,3];
			double pz = T[2,3];

			double Euler_phi = Math.Atan(ay/ax);
			double Euler_theta = Math.Atan((Math.Cos(Euler_phi)*ax+Math.Sin(Euler_phi)*ay)/az);
			double Euler_psi = Math.Atan((-Math.Sin(Euler_phi)*nx+Math.Cos(Euler_phi)*ny)/(-Math.Sin(Euler_phi)*ox+Math.Cos(Euler_phi)*oy));

			double RPY_phi = Math.Atan(ny/nx);
			double RPY_theta = Math.Atan(-nz/(Math.Cos(RPY_phi)*nx+Math.Sin(RPY_phi)*ny));
			double RPY_psi = Math.Atan((Math.Sin(RPY_phi)*ax-Math.Cos(RPY_phi)*ay)/(-Math.Sin(RPY_phi)*ox+Math.Cos(RPY_phi)*oy));

			double[] x7 = { px,py, pz, RPY_phi*180/Math.PI, RPY_theta*180/Math.PI, RPY_psi*180/Math.PI, Euler_phi*180/Math.PI,Euler_theta*180/Math.PI,Euler_psi*180/Math.PI};
			Vector<double> result = Vec.DenseOfArray(x7);

			return result;
		}
		
		public Vector<double> MZ07_Forward_Velocity_Kinematics_PoEModel(double theta1, double theta2, double theta3, double theta4, double theta5, double theta6, double omega_θ1, double omega_θ2, double omega_θ3, double omega_θ4, double omega_θ5, double omega_θ6)
		{
			var Mat = Matrix<double>.Build;	
			var Vec = Vector<double>.Build;
			Matrix<double> M = Mat.Dense(4,4);					
			M[0,0] = 1;
			M[1,0] = 0;
			M[2,0] = 0;
			M[3,0] = 0;

			M[0,1] = 0;
			M[1,1] = -1;
			M[2,1] = 0;
			M[3,1] = 0;

			M[0,2] = 0;
			M[1,2] = 0;
			M[2,2] = -1;
			M[3,2] = 0;

			M[0,3] = 50+330+45;
			M[1,3] = 0;
			M[2,3] = 345-340-73;
			M[3,3] = 1;

			double[] x1 = { 0 , 0 , 1 ,0 , 0 , 0 };
			Vector<double> S1 = Vec.DenseOfArray(x1);	
			double[] x2 = { 0 , -1,  0 , 345 ,  0 , -50 };
			Vector<double> S2 = Vec.DenseOfArray(x2);	
			double[] x3 = { 0 , -1 , 0 , 345 ,  0 , -330-50 };
			Vector<double> S3 = Vec.DenseOfArray(x3);	
			double[] x4 = { 0 ,  0 , -1 , 0 , 50+330+45 , 0 };
			Vector<double> S4 = Vec.DenseOfArray(x4);	
			double[] x5 = { 0 , -1 , 0 ,345-340 , 0 , -50-330-45 };
			Vector<double> S5 = Vec.DenseOfArray(x5);	
			double[] x6 = { 0 ,   0 ,  -1 ,	0 ,  50+330+45 ,0 };
			Vector<double> S6 = Vec.DenseOfArray(x6);	

			PoE_Algorithm PoE = new PoE_Algorithm();
			Matrix<double> AdT1 = PoE.AdT(S1[0],S1[1],S1[2],S1[3],S1[4],S1[5],theta1*Math.PI/180);
			Matrix<double> AdT2 = PoE.AdT(S2[0],S2[1],S2[2],S2[3],S2[4],S2[5],theta2*Math.PI/180);
			Matrix<double> AdT3 = PoE.AdT(S3[0],S3[1],S3[2],S3[3],S3[4],S3[5],theta3*Math.PI/180);
			Matrix<double> AdT4 = PoE.AdT(S4[0],S4[1],S4[2],S4[3],S4[4],S4[5],theta4*Math.PI/180);
			Matrix<double> AdT5 = PoE.AdT(S5[0],S5[1],S5[2],S5[3],S5[4],S5[5],theta5*Math.PI/180);
			Matrix<double> AdT6 = PoE.AdT(S6[0],S6[1],S6[2],S6[3],S6[4],S6[5],theta6*Math.PI/180);

			Matrix<double> AdT12 = AdT1*AdT2;
			Matrix<double> AdT123 = AdT12*AdT3;
			Matrix<double> AdT1234 = AdT123*AdT4;
			Matrix<double> AdT12345 = AdT1234*AdT5;

			Vector<double> Js1 = S1;
			Vector<double> Js2 = AdT1*S2;
			Vector<double> Js3 = AdT12*S3;
			Vector<double> Js4 = AdT123*S4;
			Vector<double> Js5 = AdT1234*S5;
			Vector<double> Js6 = AdT12345*S6;

			Matrix<double> Js = Mat.Dense(6,6);		
			Js[0,0] = Js1[0];
			Js[1,0] = Js1[1];
			Js[2,0] = Js1[2];
			Js[3,0] = Js1[3];
			Js[4,0] = Js1[4];
			Js[5,0] = Js1[5];

			Js[0,1] = Js2[0];
			Js[1,1] = Js2[1];
			Js[2,1] = Js2[2];
			Js[3,1] = Js2[3];
			Js[4,1] = Js2[4];
			Js[5,1] = Js2[5];

			Js[0,2] = Js3[0];
			Js[1,2] = Js3[1];
			Js[2,2] = Js3[2];
			Js[3,2] = Js3[3];
			Js[4,2] = Js3[4];
			Js[5,2] = Js3[5];

			Js[0,3] = Js4[0];
			Js[1,3] = Js4[1];
			Js[2,3] = Js4[2];
			Js[3,3] = Js4[3];
			Js[4,3] = Js4[4];
			Js[5,3] = Js4[5];

			Js[0,4] = Js5[0];
			Js[1,4] = Js5[1];
			Js[2,4] = Js5[2];
			Js[3,4] = Js5[3];
			Js[4,4] = Js5[4];
			Js[5,4] = Js5[5];

			Js[0,5] = Js6[0];
			Js[1,5] = Js6[1];
			Js[2,5] = Js6[2];
			Js[3,5] = Js6[3];
			Js[4,5] = Js6[4];
			Js[5,5] = Js6[5];

			double[] x7 = {omega_θ1,omega_θ2,omega_θ3,omega_θ4,omega_θ5,omega_θ6};
			Vector<double> omega_θ = Vec.DenseOfArray(x7);
			Vector<double> V = Js*omega_θ;

			return V;
		}
    
		public Vector<double> MZ07_Numerical_Inverse_Kinematics_PoEModel(double X_target, double Y_target, double Z_target, double Roll_target, double Pitch_target, double Yaw_target, int steps = 10, double error = 0.01, double theta0_1=10, double theta0_2=10, double theta0_3=10, double theta0_4=10, double theta0_5=10, double theta0_6=10)
		{
			//目标点计算
			Roll_target = Roll_target*Math.PI/180;
			Pitch_target = Pitch_target*Math.PI/180;
			Yaw_target = Yaw_target*Math.PI/180; 

			double r11 = Math.Cos(Pitch_target)*Math.Cos(Roll_target);
			double r12 = Math.Cos(Roll_target)*Math.Sin(Pitch_target)*Math.Sin(Yaw_target) - Math.Cos(Yaw_target)*Math.Sin(Roll_target);
			double r13 = Math.Cos(Roll_target)*Math.Cos(Yaw_target)*Math.Sin(Pitch_target) + Math.Sin(Roll_target)*Math.Sin(Yaw_target);
			double r21 = Math.Cos(Pitch_target)*Math.Sin(Roll_target);
			double r22 = Math.Cos(Roll_target)*Math.Cos(Yaw_target) + Math.Sin(Pitch_target)*Math.Sin(Roll_target)*Math.Sin(Yaw_target);
			double r23 = -Math.Cos(Roll_target)*Math.Sin(Yaw_target) + Math.Cos(Yaw_target)*Math.Sin(Pitch_target)*Math.Sin(Roll_target);
			double r31 = -Math.Sin(Pitch_target);
			double r32 = Math.Cos(Pitch_target)*Math.Sin(Yaw_target);
			double r33 = Math.Cos(Pitch_target)*Math.Cos(Yaw_target);

			PoE_Algorithm PoE = new PoE_Algorithm();
			var Mat = Matrix<double>.Build;	
			var Vec = Vector<double>.Build;

			//目标位置
			double[,] x1 = {{ r11,r12,r13,X_target },
              			    { r21,r22,r23,Y_target },
						    { r31,r32,r33,Z_target },
							{ 0,  0,  0,         1 }};
			Matrix<double> Ts_target = Mat.DenseOfArray(x1);	

			double[,] x2 = {{1,  0,  0,  50+330+45},
							{0, -1,  0,          0},
							{0,  0, -1, 345-340-73},
							{0,  0,  0,          1}};
			Matrix<double> M = Mat.DenseOfArray(x2);		

			double[,] x3 = {{ 1, 0, 0 },
              			    { 0, 1, 0 },
						    { 0, 0, 1 }};
			Matrix<double> I = Mat.DenseOfArray(x3);	

			double[] x4 = { 0 , 0 , 1 ,0 , 0 , 0 };
			Vector<double> S1 = Vec.DenseOfArray(x4);	
			double[] x5 = { 0 , -1,  0 , 345 ,  0 , -50 };
			Vector<double> S2 = Vec.DenseOfArray(x5);	
			double[] x6 = { 0 , -1 , 0 , 345 ,  0 , -330-50 };
			Vector<double> S3 = Vec.DenseOfArray(x6);	
			double[] x7 = { 0 ,  0 , -1 , 0 , 50+330+45 , 0 };
			Vector<double> S4 = Vec.DenseOfArray(x7);	
			double[] x8 = { 0 , -1 , 0 ,345-340 , 0 , -50-330-45 };
			Vector<double> S5 = Vec.DenseOfArray(x8);	
			double[] x9 = { 0 ,   0 ,  -1 ,	0 ,  50+330+45 ,0 };
			Vector<double> S6 = Vec.DenseOfArray(x9);	

			double[] x10 = { 0, 0, -1, 0, -50-330-45, 0};
			Vector<double> B1 = Vec.DenseOfArray(x10);
			double[] x11 = { 0, 1,  0, 340+73 ,0,-330-45};		
			Vector<double> B2 = Vec.DenseOfArray(x11);
			double[] x12 = { 0, 1,0,340+73, 0, -45};
			Vector<double> B3 = Vec.DenseOfArray(x12);
			double[] x13 = { 0, 0,  1,  0,  0,  0};
			Vector<double> B4 = Vec.DenseOfArray(x13);
			double[] x14 = { 0, 1 , 0, 73, 0, 0};
			Vector<double> B5 = Vec.DenseOfArray(x14);			
			double[] x15 = { 0, 0,  1,  0, 0, 0};
			Vector<double> B6 = Vec.DenseOfArray(x15);				

			List<double> errorList = new List<double>();
			List<double> theta1 = new List<double>();
			List<double> theta2 = new List<double>();
			List<double> theta3 = new List<double>();
			List<double> theta4 = new List<double>();
			List<double> theta5 = new List<double>();
			List<double> theta6 = new List<double>();

			//初始预估值位置
			double[] x16 = { theta0_1*Math.PI/180,
							 theta0_2*Math.PI/180,
						 	 theta0_3*Math.PI/180,
						 	 theta0_4*Math.PI/180,
							 theta0_5*Math.PI/180,
							 theta0_6*Math.PI/180 };
			Vector<double> theta = Vec.DenseOfArray(x16);	

			//int i = 1;
			//while(i<steps+1)
			//i = i + 1;
			for (int i = 0 ; i < steps; i++)
			{
				Matrix<double> e_Sθ1 = PoE.e_Sθ(S1[0],S1[1],S1[2],S1[3],S1[4],S1[5],theta[0]);
				Matrix<double> e_Sθ2 = PoE.e_Sθ(S2[0],S2[1],S2[2],S2[3],S2[4],S2[5],theta[1]);
				Matrix<double> e_Sθ3 = PoE.e_Sθ(S3[0],S3[1],S3[2],S3[3],S3[4],S3[5],theta[2]);
				Matrix<double> e_Sθ4 = PoE.e_Sθ(S4[0],S4[1],S4[2],S4[3],S4[4],S4[5],theta[3]);
				Matrix<double> e_Sθ5 = PoE.e_Sθ(S5[0],S5[1],S5[2],S5[3],S5[4],S5[5],theta[4]);
				Matrix<double> e_Sθ6 = PoE.e_Sθ(S6[0],S6[1],S6[2],S6[3],S6[4],S6[5],theta[5]);

				Matrix<double> Ts_b12 = e_Sθ1*e_Sθ2;
				Matrix<double> Ts_b13 = Ts_b12*e_Sθ3;
				Matrix<double> Ts_b14 = Ts_b13*e_Sθ4;
				Matrix<double> Ts_b15 = Ts_b14*e_Sθ5;
				Matrix<double> Ts_b16 = Ts_b15*e_Sθ6;
				Matrix<double> Ts_b = Ts_b16*M;

				Matrix<double> Ts_b_Inverse = Ts_b.Inverse();
				Matrix<double> Tb_s = Ts_b_Inverse;
				Matrix<double> Tb_target = Tb_s * Ts_target;

				double nx = Tb_target[0,0];
				double ny = Tb_target[1,0];
				double nz = Tb_target[2,0];
				double ox = Tb_target[0,1];
				double oy = Tb_target[1,1];
				double oz = Tb_target[2,1];
				double ax = Tb_target[0,2];
				double ay = Tb_target[1,2];
				double az = Tb_target[2,2];
				double px = Tb_target[0,3];
				double py = Tb_target[1,3];
				double pz = Tb_target[2,3];
				
				double θi = Math.Acos((nx + oy + az  - 1)/2);

				double ωi1 = (oz-ay)/(2*Math.Sin(θi));
				double ωi2 = (ax-nz)/(2*Math.Sin(θi));
				double ωi3 = (ny-ox)/(2*Math.Sin(θi));

				double[,] x17 = {{ 0, -ωi3, ωi2 },
								{ ωi3, 0, -ωi1 },
								{ -ωi2, ωi1, 0 }};			
				Matrix<double> ωi = Mat.DenseOfArray(x17);	 

				Matrix<double> Ri = I*θi+(1-Math.Cos(θi))*ωi+(θi-Math.Sin(θi))*ωi*ωi;
				Matrix<double> Ri_Inverse = Ri.Inverse();

				double[] x18 = { px ,   py ,  pz  };
				Vector<double> Pi = Vec.DenseOfArray(x18);	

				Vector<double> vi = Ri_Inverse*Pi;

				double vi1 = vi[0];
				double vi2 = vi[1];
				double vi3 = vi[2];

				double[] x19 = { ωi1*θi, ωi2*θi, ωi3*θi, vi1*θi, vi2*θi, vi3*θi };
				Vector<double> Vbi = Vec.DenseOfArray(x19);	

				Matrix<double> AdT_b1 = PoE.AdT(B1[0],B1[1],B1[2],B1[3],B1[4],B1[5],-theta[0]);
				Matrix<double> AdT_b2 = PoE.AdT(B2[0],B2[1],B2[2],B2[3],B2[4],B2[5],-theta[1]);
				Matrix<double> AdT_b3 = PoE.AdT(B3[0],B3[1],B3[2],B3[3],B3[4],B3[5],-theta[2]);
				Matrix<double> AdT_b4 = PoE.AdT(B4[0],B4[1],B4[2],B4[3],B4[4],B4[5],-theta[3]);
				Matrix<double> AdT_b5 = PoE.AdT(B5[0],B5[1],B5[2],B5[3],B5[4],B5[5],-theta[4]);
				Matrix<double> AdT_b6 = PoE.AdT(B6[0],B6[1],B6[2],B6[3],B6[4],B6[5],-theta[5]);

				Matrix<double> AdT_b65 = AdT_b6 * AdT_b5;
				Matrix<double> AdT_b654= AdT_b65 * AdT_b4;
				Matrix<double> AdT_b6543 = AdT_b654 * AdT_b3;
				Matrix<double> AdT_b65432 = AdT_b6543 * AdT_b2;

				Vector<double> Jb1 = AdT_b65432 * B1;
				Vector<double> Jb2 = AdT_b6543 * B2;
				Vector<double> Jb3 = AdT_b654 * B3;
				Vector<double> Jb4 = AdT_b65 * B4;
				Vector<double> Jb5 = AdT_b6 * B5;
				Vector<double> Jb6 = B6;

				Matrix<double> Jb = Mat.Dense(6,6);

				Jb[0,0] = Jb1[0];
				Jb[1,0] = Jb1[1];
				Jb[2,0] = Jb1[2];
				Jb[3,0] = Jb1[3];
				Jb[4,0] = Jb1[4];
				Jb[5,0] = Jb1[5];

				Jb[0,1] = Jb2[0];
				Jb[1,1] = Jb2[1];
				Jb[2,1] = Jb2[2];
				Jb[3,1] = Jb2[3];
				Jb[4,1] = Jb2[4];
				Jb[5,1] = Jb2[5];

				Jb[0,2] = Jb3[0];
				Jb[1,2] = Jb3[1];
				Jb[2,2] = Jb3[2];
				Jb[3,2] = Jb3[3];
				Jb[4,2] = Jb3[4];
				Jb[5,2] = Jb3[5];

				Jb[0,3] = Jb4[0];
				Jb[1,3] = Jb4[1];
				Jb[2,3] = Jb4[2];
				Jb[3,3] = Jb4[3];
				Jb[4,3] = Jb4[4];
				Jb[5,3] = Jb4[5];

				Jb[0,4] = Jb5[0];
				Jb[1,4] = Jb5[1];
				Jb[2,4] = Jb5[2];
				Jb[3,4] = Jb5[3];
				Jb[4,4] = Jb5[4];
				Jb[5,4] = Jb5[5];

				Jb[0,5] = Jb6[0];
				Jb[1,5] = Jb6[1];
				Jb[2,5] = Jb6[2];
				Jb[3,5] = Jb6[3];
				Jb[4,5] = Jb6[4];
				Jb[5,5] = Jb6[5];

				//求矩阵的逆
				Matrix<double> Jb_Inverse = Jb.Inverse();
				Vector<double> thetaX = theta;
				theta = theta + Jb_Inverse*Vbi;
				Vector<double> deltatheta = theta-thetaX;

				double e = Math.Sqrt(deltatheta[0]*deltatheta[0] + deltatheta[1]*deltatheta[1] + deltatheta[2]*deltatheta[2] + deltatheta[3]*deltatheta[3]+ deltatheta[4]*deltatheta[4]+ deltatheta[5]*deltatheta[5]);
				errorList.Add(e);
			
				theta1.Add(theta[0]*180/Math.PI);
				theta2.Add(theta[1]*180/Math.PI);
				theta3.Add(theta[2]*180/Math.PI);
				theta4.Add(theta[3]*180/Math.PI);
				theta5.Add(theta[4]*180/Math.PI);
				theta6.Add(theta[5]*180/Math.PI);
				if(( e < error)) 
				{
					Console.WriteLine(theta*180/Math.PI%360);		
					Console.WriteLine(i+1);				
					break;						
				}						
			}	
			double[] x20 = {theta[0]*180/Math.PI%360,theta[1]*180/Math.PI%360,theta[2]*180/Math.PI%360,theta[3]*180/Math.PI%360,theta[4]*180/Math.PI%360,theta[5]*180/Math.PI%360};	
			Vector<double> result = Vec.DenseOfArray(x20);	
			return result;
		} 
		public List<double> MZ07_Trajectory_Planning_StraightLine_TriplePolynomial(double θ1_origin,double θ2_origin,double θ3_origin,double θ4_origin,double θ5_origin,double θ6_origin,double θ1_target,double θ2_target,double θ3_target,double θ4_target,double θ5_target,double θ6_target,double T=3,double t=0.01)
		{

			PoE_Algorithm PoE = new PoE_Algorithm();
			var Mat = Matrix<double>.Build;	
			var Vec = Vector<double>.Build;

			double[,] x1 = {{1,  0,  0,  50+330+45},
							{0, -1,  0,          0},
							{0,  0, -1, 345-340-73},
							{0,  0,  0,          1}};
			Matrix<double> M = Mat.DenseOfArray(x1);		

			double[] x2 = { 0 , 0 , 1 ,0 , 0 , 0 };
			Vector<double> S1 = Vec.DenseOfArray(x2);	
			double[] x3 = { 0 , -1,  0 , 345 ,  0 , -50 };
			Vector<double> S2 = Vec.DenseOfArray(x3);	
			double[] x4 = { 0 , -1 , 0 , 345 ,  0 , -330-50 };
			Vector<double> S3 = Vec.DenseOfArray(x4);	
			double[] x5 = { 0 ,  0 , -1 , 0 , 50+330+45 , 0 };
			Vector<double> S4 = Vec.DenseOfArray(x5);	
			double[] x6 = { 0 , -1 , 0 ,345-340 , 0 , -50-330-45 };
			Vector<double> S5 = Vec.DenseOfArray(x6);	
			double[] x7 = { 0 ,   0 ,  -1 ,	0 ,  50+330+45 ,0 };
			Vector<double> S6 = Vec.DenseOfArray(x7);	

			//设定终点
			Matrix<double> e_Sθ1_target = PoE.e_Sθ(S1[0],S1[1],S1[2],S1[3],S1[4],S1[5],θ1_target*Math.PI/180);
			Matrix<double> e_Sθ2_target = PoE.e_Sθ(S2[0],S2[1],S2[2],S2[3],S2[4],S2[5],θ2_target*Math.PI/180);
			Matrix<double> e_Sθ3_target = PoE.e_Sθ(S3[0],S3[1],S3[2],S3[3],S3[4],S3[5],θ3_target*Math.PI/180);
			Matrix<double> e_Sθ4_target = PoE.e_Sθ(S4[0],S4[1],S4[2],S4[3],S4[4],S4[5],θ4_target*Math.PI/180);
			Matrix<double> e_Sθ5_target = PoE.e_Sθ(S5[0],S5[1],S5[2],S5[3],S5[4],S5[5],θ5_target*Math.PI/180);
			Matrix<double> e_Sθ6_target = PoE.e_Sθ(S6[0],S6[1],S6[2],S6[3],S6[4],S6[5],θ6_target*Math.PI/180);	


			Matrix<double> T12_target = e_Sθ1_target*e_Sθ2_target;
			Matrix<double> T13_target = T12_target*e_Sθ3_target;
			Matrix<double> T14_target = T13_target*e_Sθ4_target;
			Matrix<double> T15_target = T14_target*e_Sθ5_target;
			Matrix<double> T16_target = T15_target*e_Sθ6_target;
			Matrix<double> Ts_target = T16_target*M;

			//设定起点
			Matrix<double> e_Sθ1_origin = PoE.e_Sθ(S1[0],S1[1],S1[2],S1[3],S1[4],S1[5],θ1_origin*Math.PI/180);
			Matrix<double> e_Sθ2_origin = PoE.e_Sθ(S2[0],S2[1],S2[2],S2[3],S2[4],S2[5],θ2_origin*Math.PI/180);
			Matrix<double> e_Sθ3_origin = PoE.e_Sθ(S3[0],S3[1],S3[2],S3[3],S3[4],S3[5],θ3_origin*Math.PI/180);
			Matrix<double> e_Sθ4_origin = PoE.e_Sθ(S4[0],S4[1],S4[2],S4[3],S4[4],S4[5],θ4_origin*Math.PI/180);
			Matrix<double> e_Sθ5_origin = PoE.e_Sθ(S5[0],S5[1],S5[2],S5[3],S5[4],S5[5],θ5_origin*Math.PI/180);
			Matrix<double> e_Sθ6_origin = PoE.e_Sθ(S6[0],S6[1],S6[2],S6[3],S6[4],S6[5],θ6_origin*Math.PI/180);	


			Matrix<double> T12_origin = e_Sθ1_origin*e_Sθ2_origin;
			Matrix<double> T13_origin = T12_origin*e_Sθ3_origin;
			Matrix<double> T14_origin = T13_origin*e_Sθ4_origin;
			Matrix<double> T15_origin = T14_origin*e_Sθ5_origin;
			Matrix<double> T16_origin = T15_origin*e_Sθ6_origin;
			Matrix<double> Ts_origin = T16_origin*M;
			
			List<double> Px_Path = new List<double>();
			List<double> Py_Path = new List<double>();
			List<double> Pz_Path = new List<double>();
			List<double> vx_Path = new List<double>();
			List<double> vy_Path = new List<double>();
			List<double> vz_Path = new List<double>();
			List<double> v_Path = new List<double>();
			
			/*
			double[] px;
			double[] py;
			double[] pz;
			double[] vx;
			double[] vy;
			double[] vz;
			double[] v;
			*/

			for (double i = 0 ; i < T+t; i = i+t)
			{
				Matrix<double> Tst = Ts_origin + (3*i*i/(T*T)-2*i*i*i/(T*T*T))*(Ts_target-Ts_origin);
				Matrix<double> V_Tst = (6*i/(T*T)-6*i*i/(T*T*T))*(Ts_target-Ts_origin);

				Px_Path.Add(Tst[0,3]);
				Py_Path.Add(Tst[1,3]);
				Pz_Path.Add(Tst[2,3]);

				vx_Path.Add(V_Tst[0,3]);
				vy_Path.Add(V_Tst[1,3]);
				vz_Path.Add(V_Tst[2,3]);

				v_Path.Add(Math.Sqrt(V_Tst[0,3]*V_Tst[0,3] + V_Tst[1,3]*V_Tst[1,3] + V_Tst[2,3]*V_Tst[2,3]));		

				/*
				px[j] = Tst[0,3];
				py[j] = Tst[1,3];
				pz[j] = Tst[2,3];

				vx[j] = V_Tst[0,3];
				vy[j] = V_Tst[1,3];
				vz[j] = V_Tst[2,3];

				v[j] = Math.Sqrt(V_Tst[0,3]*V_Tst[0,3] + V_Tst[1,3]*V_Tst[1,3] + V_Tst[2,3]*V_Tst[2,3]);		
				j = j + 1;
				*/
			}

			return v_Path;
		}


		public Matrix<float> MZ07_Trajectory_Planning_StraightLine_TriplePolynomial_AngleOutput(double θ1_origin, double θ2_origin, double θ3_origin, double θ4_origin, double θ5_origin, double θ6_origin, double θ1_target, double θ2_target, double θ3_target, double θ4_target, double θ5_target, double θ6_target, double T=3, double t=0.01, int steps=100, double error=0.001)
		{

			PoE_Algorithm PoE = new PoE_Algorithm();
			var Mat = Matrix<double>.Build;	
			var MatFloat = Matrix<float>.Build;	
			var Vec = Vector<double>.Build;

			double[,] x1 = {{1,  0,  0,  50+330+45},
							{0, -1,  0,          0},
							{0,  0, -1, 345-340-73},
							{0,  0,  0,          1}};
			Matrix<double> M = Mat.DenseOfArray(x1);	

			double[,] x2 = {{ 1, 0, 0 },
							{ 0, 1, 0 },
							{ 0, 0, 1 }};
			Matrix<double> I = Mat.DenseOfArray(x2);	

			List<double> errorList = new List<double>();
			List<double> theta1 = new List<double>();
			List<double> theta2 = new List<double>();
			List<double> theta3 = new List<double>();
			List<double> theta4 = new List<double>();
			List<double> theta5 = new List<double>();
			List<double> theta6 = new List<double>();

			double[] x3 = { 0 , 0 , 1 ,0 , 0 , 0 };
			Vector<double> S1 = Vec.DenseOfArray(x3);	
			double[] x4 = { 0 , -1,  0 , 345 ,  0 , -50 };
			Vector<double> S2 = Vec.DenseOfArray(x4);	
			double[] x5 = { 0 , -1 , 0 , 345 ,  0 , -330-50 };
			Vector<double> S3 = Vec.DenseOfArray(x5);	
			double[] x6 = { 0 ,  0 , -1 , 0 , 50+330+45 , 0 };
			Vector<double> S4 = Vec.DenseOfArray(x6);	
			double[] x7 = { 0 , -1 , 0 ,345-340 , 0 , -50-330-45 };
			Vector<double> S5 = Vec.DenseOfArray(x7);	
			double[] x8 = { 0 ,   0 ,  -1 ,	0 ,  50+330+45 ,0 };
			Vector<double> S6 = Vec.DenseOfArray(x8);	

			double[] x9 = { 0, 0, -1, 0, -50-330-45, 0};
			Vector<double> B1 = Vec.DenseOfArray(x9);
			double[] x10 = { 0, 1,  0, 340+73 ,0,-330-45};		
			Vector<double> B2 = Vec.DenseOfArray(x10);
			double[] x11 = { 0, 1,0,340+73, 0, -45};
			Vector<double> B3 = Vec.DenseOfArray(x11);
			double[] x12 = { 0, 0,  1,  0,  0,  0};
			Vector<double> B4 = Vec.DenseOfArray(x12);
			double[] x13 = { 0, 1 , 0, 73, 0, 0};
			Vector<double> B5 = Vec.DenseOfArray(x13);			
			double[] x14 = { 0, 0,  1,  0, 0, 0};
			Vector<double> B6 = Vec.DenseOfArray(x14);	

			//设定终点
			double[] x15 = { θ1_target*Math.PI/180,
							 θ2_target*Math.PI/180,
							 θ3_target*Math.PI/180,
							 θ4_target*Math.PI/180,
							 θ5_target*Math.PI/180,
							 θ6_target*Math.PI/180 };
			Vector<double> θ_target = Vec.DenseOfArray(x15);	

			Matrix<double> e_Sθ1_target = PoE.e_Sθ(S1[0],S1[1],S1[2],S1[3],S1[4],S1[5],θ_target[0]);
			Matrix<double> e_Sθ2_target = PoE.e_Sθ(S2[0],S2[1],S2[2],S2[3],S2[4],S2[5],θ_target[1]);
			Matrix<double> e_Sθ3_target = PoE.e_Sθ(S3[0],S3[1],S3[2],S3[3],S3[4],S3[5],θ_target[2]);
			Matrix<double> e_Sθ4_target = PoE.e_Sθ(S4[0],S4[1],S4[2],S4[3],S4[4],S4[5],θ_target[3]);
			Matrix<double> e_Sθ5_target = PoE.e_Sθ(S5[0],S5[1],S5[2],S5[3],S5[4],S5[5],θ_target[4]);
			Matrix<double> e_Sθ6_target = PoE.e_Sθ(S6[0],S6[1],S6[2],S6[3],S6[4],S6[5],θ_target[5]);

			Matrix<double> T12_target = e_Sθ1_target*e_Sθ2_target;
			Matrix<double> T13_target = T12_target*e_Sθ3_target;
			Matrix<double> T14_target = T13_target*e_Sθ4_target;
			Matrix<double> T15_target = T14_target*e_Sθ5_target;
			Matrix<double> T16_target = T15_target*e_Sθ6_target;
			Matrix<double> Ts_target = T16_target*M;

	    	//设定起点
			double[] x16 = { θ1_origin*Math.PI/180,
							 θ2_origin*Math.PI/180,
							 θ3_origin*Math.PI/180,
							 θ4_origin*Math.PI/180,
							 θ5_origin*Math.PI/180,
							 θ6_origin*Math.PI/180 };
			Vector<double> θ_origin  = Vec.DenseOfArray(x16);	

			Matrix<double> e_Sθ1_origin = PoE.e_Sθ(S1[0],S1[1],S1[2],S1[3],S1[4],S1[5],θ_origin[0]);
			Matrix<double> e_Sθ2_origin = PoE.e_Sθ(S2[0],S2[1],S2[2],S2[3],S2[4],S2[5],θ_origin[1]);
			Matrix<double> e_Sθ3_origin = PoE.e_Sθ(S3[0],S3[1],S3[2],S3[3],S3[4],S3[5],θ_origin[2]);
			Matrix<double> e_Sθ4_origin = PoE.e_Sθ(S4[0],S4[1],S4[2],S4[3],S4[4],S4[5],θ_origin[3]);
			Matrix<double> e_Sθ5_origin = PoE.e_Sθ(S5[0],S5[1],S5[2],S5[3],S5[4],S5[5],θ_origin[4]);
			Matrix<double> e_Sθ6_origin = PoE.e_Sθ(S6[0],S6[1],S6[2],S6[3],S6[4],S6[5],θ_origin[5]);

			Matrix<double> T12_origin = e_Sθ1_origin*e_Sθ2_origin;
			Matrix<double> T13_origin = T12_origin*e_Sθ3_origin;
			Matrix<double> T14_origin = T13_origin*e_Sθ4_origin;
			Matrix<double> T15_origin = T14_origin*e_Sθ5_origin;
			Matrix<double> T16_origin = T15_origin*e_Sθ6_origin;
			Matrix<double> Ts_origin = T16_origin*M;


			List<double> Px_Path = new List<double>();
			List<double> Py_Path = new List<double>();
			List<double> Pz_Path = new List<double>();
			List<double> vx_Path = new List<double>();
			List<double> vy_Path = new List<double>();
			List<double> vz_Path = new List<double>();
			List<double> v_Path = new List<double>();

			Matrix<double> V_Tst = Mat.Dense(4,4);
			Vector<double> theta = Vec.Dense(6);

			for (double i = 0 ; i < T+t; i = i+t)
			{
				if(i == 0)
				{    
					theta = θ_origin;

					theta1.Add(theta[0]*180/Math.PI);
					theta2.Add(theta[1]*180/Math.PI);
					theta3.Add(theta[2]*180/Math.PI);
					theta4.Add(theta[3]*180/Math.PI);
					theta5.Add(theta[4]*180/Math.PI);
					theta6.Add(theta[5]*180/Math.PI);

					Px_Path.Add(Ts_origin[0,3]);
					Py_Path.Add(Ts_origin[1,3]);
					Pz_Path.Add(Ts_origin[2,3]);

					V_Tst = (6*i/(T*T)-6*i*i/(T*T*T))*(Ts_target-Ts_origin);

					vx_Path.Add(V_Tst[0,3]);
					vy_Path.Add(V_Tst[1,3]);
					vz_Path.Add(V_Tst[2,3]);

					v_Path.Add(Math.Sqrt(V_Tst[0,3]*V_Tst[0,3] + V_Tst[1,3]*V_Tst[1,3] + V_Tst[2,3]*V_Tst[2,3]));
			
					i = i + t;
				}
				Matrix<double> Tst = Ts_origin + (3*i*i/(T*T)-2*i*i*i/(T*T*T))*(Ts_target-Ts_origin);

				//逆运动学数值解
				for (int j = 0 ; j < steps; j++)
				{

					Matrix<double> e_Sθ1 = PoE.e_Sθ(S1[0],S1[1],S1[2],S1[3],S1[4],S1[5],theta[0]);
					Matrix<double> e_Sθ2 = PoE.e_Sθ(S2[0],S2[1],S2[2],S2[3],S2[4],S2[5],theta[1]);
					Matrix<double> e_Sθ3 = PoE.e_Sθ(S3[0],S3[1],S3[2],S3[3],S3[4],S3[5],theta[2]);
					Matrix<double> e_Sθ4 = PoE.e_Sθ(S4[0],S4[1],S4[2],S4[3],S4[4],S4[5],theta[3]);
					Matrix<double> e_Sθ5 = PoE.e_Sθ(S5[0],S5[1],S5[2],S5[3],S5[4],S5[5],theta[4]);
					Matrix<double> e_Sθ6 = PoE.e_Sθ(S6[0],S6[1],S6[2],S6[3],S6[4],S6[5],theta[5]);

					Matrix<double> Ts_b12 = e_Sθ1*e_Sθ2;
					Matrix<double> Ts_b13 = Ts_b12*e_Sθ3;
					Matrix<double> Ts_b14 = Ts_b13*e_Sθ4;
					Matrix<double> Ts_b15 = Ts_b14*e_Sθ5;
					Matrix<double> Ts_b16 = Ts_b15*e_Sθ6;
					Matrix<double> Ts_b = Ts_b16*M;

					Matrix<double> Ts_b_Inverse = Ts_b.Inverse();
					Matrix<double> Tb_s = Ts_b_Inverse;
					Matrix<double> Tb_target = Tb_s * Tst;		


					double nx = Tb_target[0,0];
					double ny = Tb_target[1,0];
					double nz = Tb_target[2,0];
					double ox = Tb_target[0,1];
					double oy = Tb_target[1,1];
					double oz = Tb_target[2,1];
					double ax = Tb_target[0,2];
					double ay = Tb_target[1,2];
					double az = Tb_target[2,2];
					double px = Tb_target[0,3];
					double py = Tb_target[1,3];
					double pz = Tb_target[2,3];

					double θi = Math.Acos((nx + oy + az  - 1)/2);

					double ωi1 = (oz-ay)/(2*Math.Sin(θi));
					double ωi2 = (ax-nz)/(2*Math.Sin(θi));
					double ωi3 = (ny-ox)/(2*Math.Sin(θi));

					double[,] x17 = {{ 0, -ωi3, ωi2 },
									{ ωi3, 0, -ωi1 },
									{ -ωi2, ωi1, 0 }};			
					Matrix<double> ωi = Mat.DenseOfArray(x17);	 

					Matrix<double> Ri = I*θi+(1-Math.Cos(θi))*ωi+(θi-Math.Sin(θi))*ωi*ωi;
					Matrix<double> Ri_Inverse = Ri.Inverse();

					double[] x18 = { px ,   py ,  pz  };
					Vector<double> Pi = Vec.DenseOfArray(x18);	

					Vector<double> vi = Ri_Inverse*Pi;

					double vi1 = vi[0];
					double vi2 = vi[1];
					double vi3 = vi[2];

					double[] x19 = { ωi1*θi, ωi2*θi, ωi3*θi, vi1*θi, vi2*θi, vi3*θi };
					Vector<double> Vbi = Vec.DenseOfArray(x19);	

					Matrix<double> AdT_b1 = PoE.AdT(B1[0],B1[1],B1[2],B1[3],B1[4],B1[5],-theta[0]);
					Matrix<double> AdT_b2 = PoE.AdT(B2[0],B2[1],B2[2],B2[3],B2[4],B2[5],-theta[1]);
					Matrix<double> AdT_b3 = PoE.AdT(B3[0],B3[1],B3[2],B3[3],B3[4],B3[5],-theta[2]);
					Matrix<double> AdT_b4 = PoE.AdT(B4[0],B4[1],B4[2],B4[3],B4[4],B4[5],-theta[3]);
					Matrix<double> AdT_b5 = PoE.AdT(B5[0],B5[1],B5[2],B5[3],B5[4],B5[5],-theta[4]);
					Matrix<double> AdT_b6 = PoE.AdT(B6[0],B6[1],B6[2],B6[3],B6[4],B6[5],-theta[5]);

					Matrix<double> AdT_b65 = AdT_b6 * AdT_b5;
					Matrix<double> AdT_b654= AdT_b65 * AdT_b4;
					Matrix<double> AdT_b6543 = AdT_b654 * AdT_b3;
					Matrix<double> AdT_b65432 = AdT_b6543 * AdT_b2;

					Vector<double> Jb1 = AdT_b65432 * B1;
					Vector<double> Jb2 = AdT_b6543 * B2;
					Vector<double> Jb3 = AdT_b654 * B3;
					Vector<double> Jb4 = AdT_b65 * B4;
					Vector<double> Jb5 = AdT_b6 * B5;
					Vector<double> Jb6 = B6;

					Matrix<double> Jb = Mat.Dense(6,6);

					Jb[0,0] = Jb1[0];
					Jb[1,0] = Jb1[1];
					Jb[2,0] = Jb1[2];
					Jb[3,0] = Jb1[3];
					Jb[4,0] = Jb1[4];
					Jb[5,0] = Jb1[5];

					Jb[0,1] = Jb2[0];
					Jb[1,1] = Jb2[1];
					Jb[2,1] = Jb2[2];
					Jb[3,1] = Jb2[3];
					Jb[4,1] = Jb2[4];
					Jb[5,1] = Jb2[5];

					Jb[0,2] = Jb3[0];
					Jb[1,2] = Jb3[1];
					Jb[2,2] = Jb3[2];
					Jb[3,2] = Jb3[3];
					Jb[4,2] = Jb3[4];
					Jb[5,2] = Jb3[5];

					Jb[0,3] = Jb4[0];
					Jb[1,3] = Jb4[1];
					Jb[2,3] = Jb4[2];
					Jb[3,3] = Jb4[3];
					Jb[4,3] = Jb4[4];
					Jb[5,3] = Jb4[5];

					Jb[0,4] = Jb5[0];
					Jb[1,4] = Jb5[1];
					Jb[2,4] = Jb5[2];
					Jb[3,4] = Jb5[3];
					Jb[4,4] = Jb5[4];
					Jb[5,4] = Jb5[5];

					Jb[0,5] = Jb6[0];
					Jb[1,5] = Jb6[1];
					Jb[2,5] = Jb6[2];
					Jb[3,5] = Jb6[3];
					Jb[4,5] = Jb6[4];
					Jb[5,5] = Jb6[5];

					//求矩阵的逆
					Matrix<double> Jb_Inverse = Jb.Inverse();
					Vector<double> thetaX = theta;
					theta = theta + Jb_Inverse*Vbi;
					Vector<double> deltatheta = theta-thetaX;

					double e = Math.Sqrt(deltatheta[0]*deltatheta[0] + deltatheta[1]*deltatheta[1] + deltatheta[2]*deltatheta[2] + deltatheta[3]*deltatheta[3]+ deltatheta[4]*deltatheta[4]+ deltatheta[5]*deltatheta[5]);
					errorList.Add(e);
				
					if(e < error)
					{
						theta1.Add(theta[0]*180/Math.PI);
						theta2.Add(theta[1]*180/Math.PI);
						theta3.Add(theta[2]*180/Math.PI);
						theta4.Add(theta[3]*180/Math.PI);
						theta5.Add(theta[4]*180/Math.PI);
						theta6.Add(theta[5]*180/Math.PI);

						break;
					}	
				}
				Px_Path.Add(Tst[0,3]);
				Py_Path.Add(Tst[1,3]);
				Pz_Path.Add(Tst[2,3]);

				V_Tst = (6*i/(T*T)-6*i*i/(T*T*T))*(Ts_target-Ts_origin);

				vx_Path.Add(V_Tst[0,3]);
				vy_Path.Add(V_Tst[1,3]);
				vz_Path.Add(V_Tst[2,3]);

				v_Path.Add(Math.Sqrt(V_Tst[0,3]*V_Tst[0,3] + V_Tst[1,3]*V_Tst[1,3] + V_Tst[2,3]*V_Tst[2,3]));
			}

			Matrix<double> result = Mat.Dense(6,theta1.Count);

			for(int k=0; k<theta1.Count; k++)
			{
			 	result[0,k]= theta1[k];
				result[1,k]= theta2[k];
				result[2,k]= theta3[k];
				result[3,k]= theta4[k];
				result[4,k]= theta5[k];
				result[5,k]= theta6[k];

			}
			//double转换为float
			Vector<double> MZ07Theta1 = result.Row(0);
			Vector<double> MZ07Theta2 = result.Row(1);
			Vector<double> MZ07Theta3 = result.Row(2);
			Vector<double> MZ07Theta4 = result.Row(3);
			Vector<double> MZ07Theta5 = result.Row(4);
			Vector<double> MZ07Theta6 = result.Row(5);

			List<float> MZ07Theta1Float = new List<float>(); 
			List<float> MZ07Theta2Float = new List<float>(); 
			List<float> MZ07Theta3Float = new List<float>(); 
			List<float> MZ07Theta4Float = new List<float>(); 
			List<float> MZ07Theta5Float = new List<float>(); 
			List<float> MZ07Theta6Float = new List<float>(); 

			foreach (float value in MZ07Theta1)
			{
				MZ07Theta1Float.Add(value);
			}
			foreach (float value in MZ07Theta2)
			{
				MZ07Theta2Float.Add(value);
			}
			foreach (float value in MZ07Theta3)
			{
				MZ07Theta3Float.Add(value);
			}
			foreach (float value in MZ07Theta4)
			{
				MZ07Theta4Float.Add(value);
			}
			foreach (float value in MZ07Theta5)
			{
				MZ07Theta5Float.Add(value);
			}
			foreach (float value in MZ07Theta6)
			{
				MZ07Theta6Float.Add(value);
			}

			Matrix<float> resultFloat = MatFloat.Dense(6,MZ07Theta1.Count);

			for(int o=0; o < MZ07Theta1.Count; o++)
			{
			 	resultFloat[0,o]= MZ07Theta1Float[o];
				resultFloat[1,o]= MZ07Theta2Float[o];
				resultFloat[2,o]= MZ07Theta3Float[o];
				resultFloat[3,o]= MZ07Theta4Float[o];
				resultFloat[4,o]= MZ07Theta5Float[o];
				resultFloat[5,o]= MZ07Theta6Float[o];
			}

			return resultFloat;
		}
	}
}