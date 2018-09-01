# Emergence

This is the companion repository for the following article: *Human detection of probabilistic versus deterministic regularities in sequences*. The repository contains the scripts and functions needed to reproduce the analyses, simulations and fits presented in the paper.

In order to get the toolbox up and running, run the following commands in the MATLAB command window:

```
>> system('git clone https://github.com/maheump/Emergence.git')
>> cd Emergence
>> SETUP
```

The script ```SETUP.m``` sets up the environment required to run the different scripts. It should be re-run each time MATLAB is restarted as the toolbox do not save any change (be them paths added to MATLAB path list or figures' default properties) it does to the MATLAB default configurations.

After calling that script, you are good to go!

## Abstract of the paper

The environment is highly structured in time thus allowing the brain to generate expectations about future events. To maximize the accuracy of its predictions, and thus better adjust behavior, the brain must however correctly identify the process generating the observations it receives. In particular, it must accurately distinguish between processes for which it is possible to perfectly predict the future (deterministic processes) from those for which a certain degree of uncertainty remains (probabilistic processes). How does the brain infer and arbitrate between these two hypotheses? Using a novel sequence learning experiment combined with finger-tracking, we show that human subjects are able to detect both types of regularities when suddenly appearing out of randomness. Crucially, they were able to estimate, and to update as observations were delivered, the likelihood of each of these two hypotheses in a manner that was predicted by normative principles. This suggests that the brain may approximate Bayesian inference when comparing hypotheses that involved different learning policies such as to maximize its predictability power.

## Example inference

Here is an example sequence with the beliefs from an example subject and from the ideal observer when presented with the very same sequence.
![Example inference](https://github.com/maheump/Emergence/blob/initialdev/Finger%20tracking%20analyses/figs/F_M.gif)

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
