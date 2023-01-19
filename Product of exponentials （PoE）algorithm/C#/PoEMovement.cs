using System.Collections.Generic;
using UnityEngine;
using System;
using MathNet.Numerics.LinearAlgebra;
using System.Collections;
using PoE_Kinematics;
using System.Threading;  

public class PoEMovement : MonoBehaviour
{
    private List<ArticulationBody> m_RobotJoints;
    public List<float> jointTargets;
    public List<float> jointOrigins;
    public float T;
    public float interval;
    public int steps;
    public float error;
    //public int launch;

    void Start()
    {
        m_RobotJoints = new List<ArticulationBody>();

        jointOrigins = new List<float>();
        jointTargets = new List<float>();

        // 获取机器人的所有旋转关节
        foreach (ArticulationBody joint in gameObject.GetComponentsInChildren<
        ArticulationBody>())
        {
            if (joint.jointType == ArticulationJointType.RevoluteJoint)
            {
                m_RobotJoints.Add(joint);
                jointOrigins.Add(0.0f);
                jointTargets.Add(0.0f);
            }
        }
    }
    //void Update()
   // {
    //public void JointsMoveToTargets()
        //{
            // 将每个关节的target设置为jointTargets中的值
            //for (int i = 0; i < m_RobotJoints.Count; i++)
            //{

            //}
       // }
   // }
    // 开放给外部用于单独设置jointTargets的方法
    // public void SetJointTargetByIdx(int idx, float value)
    // {
    //     jointTargets[idx] = value;
    //     var jointXDrive = m_RobotJoints[idx].xDrive;
    //     jointXDrive.target = jointTargets[idx];
    //     m_RobotJoints[idx].xDrive = jointXDrive;
    // }

    //调用轨迹插补算法

    public void JointsMoveToJoints()
    {
        PoE_Algorithm n = new PoE_Algorithm();
        //PoEMovement h = new PoEMovement();
        Matrix<float> MZ07 = n.MZ07_Trajectory_Planning_StraightLine_TriplePolynomial_AngleOutput(jointOrigins[0],jointOrigins[1],jointOrigins[2],jointOrigins[3],jointOrigins[4],jointOrigins[5],jointTargets[0],jointTargets[1],jointTargets[2],jointTargets[3],jointTargets[4],jointTargets[5],T, interval, steps, error);
        Vector<float> MZ07Theta1 = MZ07.Row(0);
        
        for(int k = 0; k < MZ07Theta1.Count; k++)
        {
            // if(k<MZ07Theta1.Count)
            // {
                for(int j = 0; j < m_RobotJoints.Count; j++)
                {
                    jointTargets[j] = MZ07[j,k];
                    var jointXDrive = m_RobotJoints[j].xDrive;
                    jointXDrive.target = jointTargets[j];
                    m_RobotJoints[j].xDrive = jointXDrive;
                }
                //Thread.Sleep(16);
            //}
            
            
        }
    }
    // void Update()
    // {
    //     // void Litchigo()
    //     // {
    //         if(launch==1)
    //         {
    //                 int k = 0;
    //                 JointsMoveToJoints(k++);
    //                 Thread.Sleep(33);
    //             }
    //         //}
    // }


}

