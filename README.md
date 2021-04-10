# Simple-Inverse-Dynamic
The repository contains a script for solving the inverse problem of dynamics using the RNEA (Recursive Newton-Euler Algorithm). The script contains a description of the model, as well as a method for solving it, implemented in C#.
The project uses the MathNet.Numerics library. Kinematic and dynamic parameters of the kuka_kr6r700 manipulator are introduced in the default constructor.

You can read more about the method:
- http://royfeatherstone.org/spatial/v2/
- the book Rigid Body Dynamics Algorithms (RBDA)
- Springer Handbook of Robotics

![Algorithm](https://github.com/zik-888/Simple-Inverse-Dynamic/blob/main/Algorithm.PNG "Coordinate-free recursive Newton–
Euler algorithm (RNEA) for inverse dynamics")
