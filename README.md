# Emergence

This is the companion repository for the following article: *Rational arbitration between statistics and rules in human sequence learning*. The repository contains the scripts and functions needed to reproduce the analyses, simulations and fits presented in the paper.

## Abstract of the paper

Detecting and learning temporal regularities is essential to accurately predict the future. Past research indicates that humans are sensitive to two types of sequential regularities: deterministic rules, which afford sure predictions, and statistical biases, which govern the probabilities of individual items and their transitions. Yet, how does the human brain arbitrate between those two systems? We used finger tracking to continuously monitor the online buildup of evidence, confidence, false alarms and changes-of-mind during sequence learning. All these aspects of behaviour conformed tightly to a hierarchical Bayesian inference model with distinct hypothesis spaces for statistics vs. rules, yet linked by a single probabilistic currency. Alternative models based either on a single statistical mechanism or on two non-commensurable systems were rejected. Our results indicate that a unified Bayesian inference mechanism, capable of operating over several distinct hypothesis spaces, underlies the human capability to learn both statistics and rules.

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

These are the toolboxes and functions the toolbox depends upon. They are automatically downloaded and installed in the Dependencies folder of the toolbox.

* [Custom functions](https://github.com/maheump/matlab/)

## Related toolbox

We developed another related toolbox that features Bayesian (near-)ideal observers estimating different statistics from an input sequence. It can, for instance, deal with sequences entailing more than a single change point.

* [Minimal transition probability model](https://github.com/florentmeyniel/MinimalTransitionProbsModel)
