# conditionally_stable_circuits
This is the supporting notebook and code to the paper titled **A feedback loop of conditionally stable circuits drives the cell cycle from checkpoint to checkpoint**, submitted to Scientific Reports and available on bioRxiv under https://www.biorxiv.org/content/10.1101/654863v1

# Supporting Information - Jupyter Notebook

## A feedback loop of conditionally stable circuits drives the cell cycle from checkpoint to checkpoint
Dávid Deritei<sup>1,2</sup>, Jordan Rozum <sup>1</sup>, Erzsébet Ravasz Regan<sup>3</sup>, Réka Albert <sup>1,*</sup>

<sup>1</sup> - Department of Physics, Pennsylvania State University, University Park, PA, United States of America<br>
<sup>2</sup> - Department of Network and Data Science, Central European University, Budapest, Hungary<br>
<sup>3</sup> - Biochemistry and Molecular Biology, The College of Wooster, Wooster, OH, United States of America <br>

<sup>*</sup> - Corresponding author<br>
<br>

### Abstract

We perform logic-based network analysis on a model of the mammalian cell cycle. This model is composed of a Restriction Switch driving cell cycle commitment and a Phase Switch driving mitotic entry and exit. By generalizing the concept of stable motif, i.e., a self-sustaining positive feedback loop that maintains an associated state, we introduce the concept of conditionally stable motif, the stability of which is contingent on external conditions. We show that the stable motifs of the Phase Switch are contingent on the state of three nodes through which it receives input from the rest of the network. Biologically, these conditions correspond to cell cycle checkpoints. Holding these nodes locked (akin to a checkpoint-free cell) transforms the Phase Switch into an autonomous oscillator that robustly toggles through the cell cycle phases G1, G2 and mitosis. The conditionally stable motifs of the Phase Switch Oscillator are organized into an ordered sequence, such that they serially stabilize each other but also cause their own destabilization. Along the way they channel the dynamics of the module onto a narrow path in state space, lending robustness to the oscillation. Self-destabilizing conditionally stable motifs suggest a general negative feedback mechanism leading to sustained oscillations.

### The notebook

The **Supplementary_notebook.ipynb** Jupyter Notebook runs Python 2.7 and contains the steps and procedures to reproduce the key findings of our paper (under review, also availale on bioRxiv: https://www.biorxiv.org/content/10.1101/654863v1).
The notebook doesn't reproduce *all* the results published in the paper, however it provides the tools for anyone interested in doing similar studies, mostly in form of functions.
We will keep updating it with new implementations as soon as we have them.

The Boolean models and supporting files are included in the repository.

The necessary python packages to run the notebook:<br>
numpy, matplotlib, seaborn (optional)<br>

BooleanNet - https://github.com/ialbert/booleannet <br>
NetworkX - https://networkx.github.io/

The simulations also use the BooleanDOI (https://github.com/yanggangthu/BooleanDOI) package but in this case a working version is included locally.

Every other function used in the notebook is defined and imported from the cool_bool_tools.py file. <br>

Comments and suggestions are welcome!
