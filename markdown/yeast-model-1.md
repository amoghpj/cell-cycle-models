We have now put together a complete cell-cycle machine: the cell
starts out in G1 phase (high Cdh1), and as the mass of the cell increases,
the Cdc20 is activated. As a consequence, CycB rises rapidly, the cell is
drawn to S/G2/M, leading to an increase in Cdh1. The cell finally divides,
dropping the mass back to its initial value, and the G1 configuration is 
readopted.

In the budding yeast, a few other molecular players contribute to this story.
The first, is the CycB/Cdk1 inhibitor/binding partner, denoted the *C*yclin
dependent-*K*inase *I*nhibitor, or CKI. This molecule initially binds, and inactivates
CycB/Cdk1 activity, but this inhibition decreases as the mass of the cell increases.
Simultaneously, the cell "commits" to S-phase, known as the START transition; this
is signaled by the activity of a "Starter Kinases" denoted SK, which are A-type cyclins,
expressed in G1. (The expression of these cyclins, Cln1-2 in budding yeast,
 on specific cell-cycle dependent transcription factors, SBF in budding yeast.)
Please see Tyson and Novak, 2001 for a detailed discussion.

Now that we have a complete yeast cell cycle model, we can start exploring
results from genetic perturbation experiments. If the components that compose
our model are indeed reflective of the underlysing biology, then we should be 
able to "simulate" gene deletion phenotypes.

Below is one such example: Consider a *cln1 cln2* strain, lacking the Starter Kinase.
How do we expect such a mutant to behave? Examining the figure above, a SK-deletion
would mean that CKI is never phosphorylated/degraded, which will mean that the cell
will remain in G1 phase. In the model such a deletion would be represented by setting
the parameter $k_{13}=0$. 

Use the menu below to simulate an "SK-deletion", or a *cln1 cln2* strain.
