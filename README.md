# Numerical simulation of a hydraulic valve actuator

  Goal of this project is to simulate the movement of the piston and get preassure values in the hidraulic fluid of a hidraulic valve actuator.

Code for numerical simulation is written in C++ using the standard library and [boost](https://www.boost.org/) library (for [odeint](http://headmyshoulder.github.io/odeint-v2/)).
After the simulation python code is inovked for a visual representation of data.


Schematic of the valve in question:

<img src=images/Valve.png height="350" />

Governing equations for hidraulic fluid preassure are gives as follows:

<img src=images/jedTlak.png height="50" width="200" />

Which can be written in the discrete form (dxdt[0] in [code](main.cpp)) as:

<img src=images/discreteTlak.png height="50" />

Also equations for acceleration can b written in the discrete form (dxdt[2] in [code](main.cpp)) as:

<img src=images/discreteAcceleration.png height="50" />

Input values for the system are given as flow of the hidraulic fluid coming into the cylinder:

<img src=images/Input.png height="150" width="650"/>

Expexted output of the program should be the console output reporting the number of integration steps, and a matplotlib window with graphs.

<img src=images/output.png height="300" />
