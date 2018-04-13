# Emergence

This is the companion repository for the following article: *Human detection of probabilistic versus deterministic regularities in sequences: same inference mechanism but distinct hypothesis spaces*.

The repository contains the codes to reproduce the analyses, simulations and fits presented in the paper.

The code is written for MATLAB, it has been developed under MATLAB 2018a.

## Abstract of the paper

## Example inference

Here is an example sequence with the beliefs from an example subject and from the ideal observer when presented with the very same sequence.
![Example inference](https://user-images.githubusercontent.com/5986212/38567248-fbe65660-3ce5-11e8-957e-0116b8eafc00.gif)

## Organization of the repository

* **Ideal observer**: Functions implementing the full ideal observer of the task as well. "Toy examples" simulation scripts are also provided to easily run the different models and inspect their behavior based on visual displays.
* **Finger tracking analyses**:
* **Functions**: A set of functions used by many different scripts in the repository. It contains mostly

## Dependencies

* [Custom functions for sequence processing](https://github.com/maheump/matlab/tree/master/sequences)
* [VBA toolbox](http://mbb-team.github.io/VBA-toolbox/)
* [avconv](https://libav.org/avconv.html)

## Related toolboxes

* [Minimal transition probability model](https://github.com/florentmeyniel/MinimalTransitionProbsModel)
