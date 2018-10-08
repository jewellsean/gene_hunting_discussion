# Discussion of "Gene hunting with hidden Markov model knockoffs" by Jewell and Witten 

All figures/simulations in the discussion of Sesia et al (2018) by Jewell and Witten can be reproduced with the script ``create_jewell_witten_discussion_figures.sh``. 

In particular, the script reproduces 

 - Figure 1 

 and a (small) subset of the experiments required to reproduce 

 - Figure 2 

 To exactly reproduce Figure 2, change the ``SIGNAL`` parameter in ```run_experiments.sh``` and ``number_simulations`` in ```configs/config_desktop.yml```. We recommend running larger scale experiments in a distributed way. 