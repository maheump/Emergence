# Emergence

This is the companion repository for the following article: *Human detection of probabilistic versus deterministic regularities in sequences*.

The repository contains the codes to reproduce the analyses, simulations and fits presented in the paper.

The code is written for MATLAB, it has been developed under MATLAB 2018a.

## Abstract of the paper

The environment is highly structured in time thus allowing the brain to generate expectations about future events. To maximize the accuracy of its predictions, and thus better adjust behavior, the brain must however correctly identify the process generating the observations it receives. In particular, it must accurately distinguish between processes for which it is possible to perfectly predict the future (deterministic processes) from those for which a certain degree of uncertainty remains (probabilistic processes). How does the brain infer and arbitrate between these two hypotheses? Using a novel sequence learning experiment combined with finger-tracking, we show that human subjects are able to detect both types of regularities when suddenly appearing out of randomness. Crucially, they were able to estimate, and to update as observations were delivered, the likelihood of each of these two hypotheses in a manner that was predicted by normative principles. This suggests that the brain may approximate Bayesian inference when comparing hypotheses that involved different learning policies such as to maximize its predictability power.

## Example inference

Here is an example sequence with the beliefs from an example subject and from the ideal observer when presented with the very same sequence.
![Example inference](https://github.com/maheump/Emergence/blob/master/Finger%20tracking%20analyses/figs/F_M.gif)

## Organization of the repository

* **Ideal observer**: Functions implementing the full ideal observer of the task as well as partial observers learning one type of regularity (e.g. Bernoulli, ...). "Toy examples" simulation scripts are also provided to easily run the different models and inspect their behavior based on visual displays.
* **Finger tracking analyses**: Series of scripts, that can be run as a pipeline, that reproduce the analyses and figures of the paper.
* **Functions**: A set of functions used by many different scripts in the repository. It contains mostly
* **Design**: A set of scripts used to explain the experimental design.

## Dependencies

* [Custom functions for sequence processing](https://github.com/maheump/matlab/tree/master/sequences)
* [VBA toolbox](http://mbb-team.github.io/VBA-toolbox/)
* [cbrewer2](https://github.com/scottclowe/cbrewer2)
* [avconv](https://libav.org/avconv.html)
* [Perceptually uniform colormaps](https://fr.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps)

## Related toolboxes

* [Minimal transition probability model](https://github.com/florentmeyniel/MinimalTransitionProbsModel)
