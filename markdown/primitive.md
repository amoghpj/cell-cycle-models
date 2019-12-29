We are nearly there! To make our model a little more 
mechanistic, we will introduce two more variables

1. Cdc20$_A$ will be the "active" form of Cdc20. In the model,
   this will interact with Cdh1, instead of the Cdc20$_T$
2. We introduce a hypothetical intermediary enzyme IEP, which
   is required to introduce the delay seen in the rise of 
   Cdc20$_T$.

Finally, we model the increase in mass as a logistic function
in order to smoothly vary the cell mass over a cell's life time.
Putting everything together, a primitive model of the cell cycle
is in place! This model demonstrates oscillations, under the condition
that the mass divides every time the cell completes mitosis, represented
here as when CycB drops below 0.1.
