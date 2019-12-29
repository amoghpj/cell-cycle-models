In the previous section, we had introduced an artificial 
variable *A* that controls the exit from mitosis. This 
variable acted as an activator of the Cdh1/APC activity.
The Cdc14 phosphatase plays this role, which is activated
indirectly by Cdc20/APC. Tyson and Novak simplify this
mechanism, by directly considering the Cdc20 activity.
The following equation models the rise in Cdc20/APC
activity in the S/G2/M phase by introducing a Hill-like
dependence on CycB activity.

$\frac{d[\text{Cdc20}_T]}{dt} = k_5' + k_5'' \frac{([\text{CycB}]m J_5)^n}{1 + ([\text{CycB}]m J_5)^n} - k_6[\text{Cdc20}_T]$

Let's step back and take stock of what we have seen so far: The CycB
nullcline is only a function of Cdh1. The Cdh1 nullcline, however, is
a function of the *mass* and the activity of Cdc20 (called *A* in the previous section). 
We have now introduced a Cdc20 equation, which ties everythin together.
Now, instead of plotting three nullclines, Tyson and Novak use the Goldbeter-Koshland
expression to express Cdh1 as a function of CycB. So we have two expressions again,
The nullclines of CycB and Cdc20$_T$, and the only free parameter
is *mass*. (Since it is  computationally expensive to recompute these nullclines
for all mass values, a low and a high value of mass have been used to 
precompute the curves below.)

As you increase the mass from 0.4 to 1.0, notice the CycB nullcline
retract from the Cdh1 nullcline, towards the right. The G1 steady
state disappears, and the cell zooms to S/G2/M.
