## The CycB/Cdk1-Cdh1/APC antagonism
We'll start at the end, the exit from mitosis. 

There are two players that coordinate the end of cell division, and the
entry into G1, namely the CycB-bound Cdk1 [kinase](https://www.uniprot.org/keywords/KW-0418), 
and the Cdh1 bound APC [ubiquitin ligase](https://www.sciencedirect.com/topics/neuroscience/anaphase-promoting-complex)
The rise in activity of the CycB/Cdk1 kinase pushes the 
cell into mitosis, spcifically into the S/G2/M phases. CycB/Cdk1
phosphorylate Cdh1 and prevent its binding to APC. In opposition,
Cdh1/APC target CycB/Cdk1 for degradation. The switch from S/G2/M to
G1 phase happens when CycB/Cdk1 is degraded, and Cdh1/APC activity is high.


These opposing forces can be modeled using a pair of non-linear ODEs, with kinetic
parameters: the *k*s are rate constants, and the *J*s are the Michaelis constants.

$\frac{d[\text{CycB}]}{dt} = k_1 - (k_2' + k_2'' [\text{Cdh1}])[\text{CycB}]$

$\frac{d[\text{Cdh1}]}{dt} = \frac{(k_3' + k_3'' A)(1- [\text{Cdh1}])}{J_3 + 1 - [\text{Cdh1}]} - \frac{k_4 m [\text{CycB}] [\text{Cdh1}]}{J_4 + [\text{Cdh1}]}$

Notice the *m* in the Cdh1 equation: it denotes the *mass* of the cell. 
The numerical solution of the above equations, the timecourses of Cdh1 (solid) and CycB (dashed),
 are plotted below. Notice the logarithmic scale of CycB activity.
The default values of these sliders (high Cdh1, low CycB, low mass) indicate that
the cell has just divided. Try increasing the mass to, and beyond, a value of 0.53.
To get a sense of how the mass affects the *dynamics* of these two opposing 
molecular factors, use the sliders to set the initial conditions of Cdh1, and CycB, and
vary the mass to see where these requlators end up at the end of the time course.
