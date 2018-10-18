# Emergence

This is the companion repository for the following article: *Human detection of probabilistic versus deterministic regularities in sequences*. The repository contains the scripts and functions needed to reproduce the analyses, simulations and fits presented in the paper.

## Abstract of the paper

The environment is highly structured in time thereby allowing the brain to generate expectations about future events. To maximize its prediction accuracy, the brain must, however, identify if the hidden process generating the observations it receives is deterministic – observations are entirely predictable in that case (e.g. traffic lights) – or probabilistic – for which a certain degree of uncertainty remains (e.g. weather). How does the brain learn temporal regularities while arbitrating between these two possible generative hypotheses? Using a novel sequence learning experiment, combined with finger tracking, we show that human subjects detect suddenly appearing probabilistic and deterministic regularities with very different learning dynamics. They are also able to estimate, and keep updated, the likelihood of these two generative hypotheses given the received observations in a manner that follows normative principles. This suggests that the brain may approximate Bayesian inference when comparing hypotheses that involve different learning policies such as to maximize its prediction accuracy.

## Example inference

Here is an example sequence with the beliefs from an example subject and from the ideal observer when presented with the very same sequence.

<p align="center">
  <img src="Finger%20tracking%20analyses/figs/F_M.gif" width="500" align="middle">
</p>

## Installation

In order to get the toolbox up and running, run the following commands in the MATLAB command window:

```
>> !git clone https://github.com/maheump/Emergence.git
>> cd Emergence
>> SETUP
```

The script ```SETUP.m``` sets up the environment required to run the different scripts. It should be re-run each time MATLAB is restarted as the toolbox do not save any change (be them paths added to MATLAB path list or figures' default properties) it does to the MATLAB default configurations.

After calling that script, you are good to go!

## Organization of the repository

* **Design**: A set of scripts used to explain the experimental design.
* **Dependencies**: Contains the general purpose MATLAB functions the toolbox depends upon.
* **Finger tracking analyses**: Series of scripts, that can be run as a pipeline, that reproduce the analyses and figures of the paper.
* **Functions**: A set of functions used by many different scripts in the repository. It contains mostly
* **Ideal observer**: Functions implementing the full ideal observer of the task as well as partial observers learning one type of regularity (e.g. Bernoulli, ...). "Toy examples" simulation scripts are also provided to easily run the different models and inspect their behavior based on visual displays.
* **Stimulation**: Scripts and functions used to run the behavioral experiment.

## Compatibility issues

The code is written for MATLAB. Note that it has been developed under MATLAB 2018a (v9.5) with the *Image Processing* (v10.3), *Signal Processing* (v8.1) and *Statistics and Machine Learning* (v11.4) toolboxes.

Compatibility with earlier versions of MATLAB has not been tested and is therefore not granted. In particular at line 151 of the function ```Emergence_IO_Tree```, instead of using the time-consuming```repmat``` function, we use a MATLAB fairly recently developed matrix logical computation.

## Dependencies

These are the toolboxes and functions the toolbox depends upon. They are automatically downloaded and installed (if not already available in the case of VBA) in the Dependencies folder of the toolbox.

* [VBA toolbox](http://mbb-team.github.io/VBA-toolbox/)
* [cbrewer2](https://github.com/maheump/matlab/)
* [Custom functions](https://github.com/maheump/matlab/)

## Related toolbox

We developed another related toolbox that features Bayesian (near-)ideal observers estimating different statistics from an input sequence. It can, for instance, deal with sequences entailing more than a single change point.

* [Minimal transition probability model](https://github.com/florentmeyniel/MinimalTransitionProbsModel)
