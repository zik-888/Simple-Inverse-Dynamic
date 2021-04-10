using System;
using System.Linq;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Dynamic
{
    public class Model
    {
        private int nb;
        private JType[] typesJoint;
        private Matrix<double>[] I;
        private Matrix<double>[] XTree;


        public Model()
        {
            nb = 7;

            typesJoint = new JType[]
            {
                JType.Rz,
                JType.Ry,
                JType.Ry,
                JType.Rx,
                JType.Ry,
                JType.Rx,
                JType.Rz
            };


            I = new Matrix<double>[7];

            I[0] = DenseMatrix.OfArray(new double[,]
            {
                {0.0334, 0, 0, 0, 0.1852, -0.0890},
                {0, 0.0284, -0.0013, -0.1852, 0, 0},
                {0, -0.0013, 0.0091, 0.0890, 0, 0},
                {0, -0.1852, 0.0890, 2.7400, 0, 0},
                {0.1852, 0, 0, 0, 2.7400, 0},
                {-0.0890, 0, 0, 0, 0, 2.7400}
            });

            I[1] = DenseMatrix.OfArray(new double[,]
            {
                {0.034, 0, 0, 0, -0.0883, -0.1880},
                {0, 0.0091, 0.0013, 0.0883, 0, 0},
                {0, 0.0013, 0.0290, 0.1880, 0, 0},
                {0, 0.0883, 0.1880, 2.7400, 0, 0},
                {-0.0883, 0, 0, 0, 2.7400, 0},
                {-0.1880, 0, 0, 0, 0, 2.7400}
            });

            I[2] = DenseMatrix.OfArray(new double[,]
            {
                {0.0125, -0.0060, 2.2413e-04, 0, 0.0755, 0.0753},
                {-0.0060, 0.0175, 2.0925e-04, -0.0755, 0, -0.1118},
                {2.2413e-04, 2.0925e-04, 0.0158, -0.0753, 0.1118, 0},
                {0, -0.0755, -0.0753, 2.3800, 0, 0},
                {0.0755, 0, 0.1118, 0, 2.3800, 0},
                {0.0753, -0.1118, 0, 0, 0, 2.3800}
            });

            I[3] = DenseMatrix.OfArray(new double[,]
            {
                {0.0131, 0.0064, 3.2280e-04, 0, -0.0759, 0.0802},
                {0.0064, 0.0138, -1.9827e-04, 0.0759, 0, 0.0858},
                {3.2280e-04, -1.9827e-04, 0.0160, -0.0802, -0.0858, 0},
                {0, 0.0759, -0.0802, 2.3800, 0, 0},
                {-0.0759, 0, -0.0858, 0, 2.3800, 0},
                {0.0802, 0.0858, 0, 0, 0, 2.3800}
            });

            I[4] = DenseMatrix.OfArray(new double[,]
            {
                {0.0703, 6.5028e-07, -1.0513e-05, 0, 0.2854, 0.1673},
                {6.5028e-07, 0.0586, 0.0097, -0.2854, 0, 0},
                {-1.0513e-05, 0.0097, 0.0147, -0.1673, 0, 0},
                {0, -0.2854, -0.1673, 2.7400, 0, 0},
                {0.2854, 0, 0, 0, 2.7400, 0},
                {0.1673, 0, 0, 0, 0, 2.7400}
            });

            I[5] = DenseMatrix.OfArray(new double[,]
            {
                {0.0033, -0.0012, -2.1222e-04, 0, -0.0165, 0.0141},
                {-0.0012, 0.0083, -1.9657e-05, 0.0165, 0, -0.0791},
                {-2.1222e-04, -1.9657e-05, 0.0098, -0.0141, 0.0791, 0},
                {0, 0.0165, -0.0141, 1.5500, 0, 0},
                {-0.0165, 0, 0.0791, 0, 1.5500, 0},
                {0.0141, -0.0791, 0, 0, 0, 1.5500}
            });

            I[6] = DenseMatrix.OfArray(new double[,]
            {
                {0.0032, -1.8634e-04, -3.4540e-04, 0, -0.0351, 0.0058},
                {-1.8634e-04, 0.0032, -4.7258e-04, 0.0351, 0, -0.0059},
                {-3.4540e-04, -4.7258e-04, 8.4561e-04, -0.0058, 0.0059, 0},
                {0, 0.0351, -0.0058, 0.5400, 0, 0},
                {-0.0351, 0, 0.0059, 0, 0.5400, 0},
                {0.0058, -0.0059, 0, 0, 0, 0.5400}
            });

            XTree = new Matrix<double>[7];

            XTree[0] = DenseMatrix.OfArray(new double[,]
            {
                {1, 0, 0, 0, 0, 0},
                {0, 1, 0, 0, 0, 0},
                {0, 0, 1, 0, 0, 0},
                {0, 0.4, 0, 1, 0, 0},
                {-0.4, 0, 0, 0, 1, 0},
                {0, 0, 0, 0, 0, 1}
            });

            XTree[1] = DenseMatrix.OfArray(new double[,]
            {
                {1, 0, 0, 0, 0, 0},
                {0, 1, 0, 0, 0, 0},
                {0, 0, 1, 0, 0, 0},
                {0, 0, 0, 1, 0, 0},
                {0, 0, 0.025, 0, 1, 0},
                {0, -0.025, 0, 0, 0, 1}
            });

            XTree[2] = DenseMatrix.OfArray(new double[,]
            {
                {1, 0, 0, 0, 0, 0},
                {0, 1, 0, 0, 0, 0},
                {0, 0, 1, 0, 0, 0},
                {0, 0, 0, 1, 0, 0},
                {0, 0, 0.315, 0, 1, 0},
                {0, -0.315, 0, 0, 0, 1}
            });

            XTree[3] = DenseMatrix.OfArray(new double[,]
            {
                {1, 0, 0, 0, 0, 0},
                {0, 1, 0, 0, 0, 0},
                {0, 0, 1, 0, 0, 0},
                {0, 0.035, 0, 1, 0, 0},
                {-0.035, 0, 0, 0, 1, 0},
                {0, 0, 0, 0, 0, 1}
            });

            XTree[4] = DenseMatrix.OfArray(new double[,]
            {
                {1, 0, 0, 0, 0, 0},
                {0, 1, 0, 0, 0, 0},
                {0, 0, 1, 0, 0, 0},
                {0, 0, 0, 1, 0, 0},
                {0, 0, 0.365, 0, 1, 0},
                {0, -0.365, 0, 0, 0, 1}
            });

            XTree[5] = DenseMatrix.OfArray(new double[,]
            {
                {1, 0, 0, 0, 0, 0},
                {0, 1, 0, 0, 0, 0},
                {0, 0, 1, 0, 0, 0},
                {0, 0, 0, 1, 0, 0},
                {0, 0, 0.08, 0, 1, 0},
                {0, -0.08, 0, 0, 0, 1}
            });

            XTree[6] = DenseMatrix.OfArray(new double[,]
            {
                {1, 0, 0, 0, 0, 0},
                {0, 1, 0, 0, 0, 0},
                {0, 0, 1, 0, 0, 0},
                {0, 0, 0, 1, 0, 0},
                {0, 0, 0, 0, 1, 0},
                {0, 0, 0, 0, 0, 1}
            });
        }


        public int NB => nb;

        public Matrix<double> GetI(int i)
        {
            return I[i];
        }

        public JType GetJType(int i)
        {
            return typesJoint[i];
        }

        public Matrix<double> GetXTree(int i)
        {
            return XTree[i];
        }

    }

    public enum JType
    {
        Rx,
        Ry,
        Rz
    }


    public class ID
    {
        public static double[] GetMethod(Model model, double[] q, double[] dq, double[] ddq)
        {
            var aGrav = GetGravity();

            var s = new Vector<double>[model.NB];
            var xup = new Matrix<double>[model.NB];

            var v = new Vector<double>[model.NB];
            var a = new Vector<double>[model.NB];
            var f = new Vector<double>[model.NB];

            Vector<double> tau = new DenseVector(model.NB);

            for (var i = 0; i < model.NB; i++)
            {
                var (XJ, S) = JCalc(model.GetJType(i), q[i]);

                s[i] = S;

                var vJ = s[i] * dq[i];

                xup[i] = XJ * model.GetXTree(i);


                if (i == 0)
                {
                    v[i] = vJ;
                    a[i] = xup[i] * -aGrav + s[i] * ddq[i];
                }
                else
                {
                    v[i] = xup[i] * v[i - 1] + vJ;
                    a[i] = xup[i] * a[i - 1] + s[i] * ddq[i] + Crm(v[i]) * vJ;
                }

                f[i] = model.GetI(i) * a[i] + Crf(v[i]) * model.GetI(i) * v[i];
            }

            for (var i = model.NB - 1; i >= 0; i--)
            {
                tau[i] = s[i] * f[i];

                if (i != 0) 
                    f[i - 1] = f[i - 1] + xup[i].Transpose() * f[i];
            }

            foreach (var VARIABLE in f)
            {
                Console.WriteLine(VARIABLE);
            }

            return tau.ToArray();
        }

        private static (Matrix<double>, Vector<double>) JCalc(JType type, double q)
        {
            var s = new DenseVector(6);

            var xj = RotX(q);

            switch (type)
            {
                case JType.Rx:
                    xj = RotX(q);
                    s[0] = 1f;
                    break;
                case JType.Ry:
                    xj = RotY(q);
                    s[1] = 1f;
                    break;
                case JType.Rz:
                    xj = RotZ(q);
                    s[2] = 1f;
                    break;
            }

            return (xj, s);
        }

        private static Matrix<double> RotX(double theta)
        {
            var c = (double) Math.Cos(theta);
            var s = (double) Math.Sin(theta);

            var m = DenseMatrix.OfArray(new double[,]
            {
                {1, 0, 0, 0, 0, 0},
                {0, c, s, 0, 0, 0},
                {0, -s, c, 0, 0, 0},
                {0, 0, 0, 1, 0, 0},
                {0, 0, 0, 0, c, s},
                {0, 0, 0, 0, -s, c}
            });

            return m;
        }

        private static Matrix<double> RotY(double theta)
        {
            var c = (double) Math.Cos(theta);
            var s = (double) Math.Sin(theta);

            var m = DenseMatrix.OfArray(new double[,]
            {
                {c, 0, -s, 0, 0, 0},
                {0, 1, 0, 0, 0, 0},
                {s, 0, c, 0, 0, 0},
                {0, 0, 0, c, 0, -s},
                {0, 0, 0, 0, 1, 0},
                {0, 0, 0, s, 0, c}
            });

            return m;
        }

        private static Matrix<double> RotZ(double theta)
        {
            var c = (double) Math.Cos(theta);
            var s = (double) Math.Sin(theta);

            var m = DenseMatrix.OfArray(new double[,]
            {
                {c, s, 0, 0, 0, 0},
                {-s, c, 0, 0, 0, 0},
                {0, 0, 1, 0, 0, 0},
                {0, 0, 0, c, s, 0},
                {0, 0, 0, -s, c, 0},
                {0, 0, 0, 0, 0, 1}
            });

            return m;
        }

        private static Vector<double> GetGravity()
        {
            return DenseVector.OfArray(new double[] {0, 0, 0, 0, 0, -9.81});
        }

        private static Matrix<double> Crm(Vector<double> v)
        {
            if (v.Count == 6)
            {
                var m = DenseMatrix.OfArray(new double[,]
                {
                    {0, -v[2], v[1], 0, 0, 0},
                    {v[2], 0, -v[0], 0, 0, 0},
                    {-v[1], v[0], 0, 0, 0, 0},
                    {0, -v[5], v[4], 0, -v[2], v[1]},
                    {v[5], 0, -v[3], v[2], 0, -v[0]},
                    {-v[4], v[3], 0, -v[1], v[0], 0}
                });

                return m;
            }
            else
            {
                var m = DenseMatrix.OfArray(new double[,]
                {
                    {0, 0, 0},
                    {v[2], 0, -v[0]},
                    {-v[1], v[0], 0}
                });

                return m;
            }
        }

        private static Matrix<double> Crf(Vector<double> v)
        {
            return -Crm(v).Transpose();
        }
    }
}