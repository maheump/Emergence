# Ideal observer

## Scripts and functions

The different observers available are:
* ```Emergence_IO_Bernoulli``` estimates the frequency of A and B.
* ```Emergence_IO_Markov``` estimates the frequency of first-order transitions (A|A, A|B, B|A, and B|B).
* ```Emergence_IO_Chain``` estimates the frequency of transitions of any order.
* ```Emergence_IO_Tree``` detects repetition of pattern of any length up to a given limit.

The full Bayesian ideal observer of the task is implemented in  ```Emergence_IO_FullIO``` and considers there might be a change point in the sequence separating a fully-stochastic part from a regular part that can be described using one of the previous observers.

Toy examples scripts are available for each of there observers.

## Full Bayesian ideal observer of the task

## Example inference

Here is the result of the inference by the full Bayesian ideal observer of the task (the function ```Emergence_IO_FullIO```).
These figures can be reproduced using the script ```Emergence_IO_ToyExampleFullIO```.

![Example inference](https://github.com/maheump/Emergence/blob/initialdev/Ideal%20observer/ToyExamples/figs/Emergence_IO_ToyExampleFullIO.jpeg)

Let's now take a look at the different subplots

The first group of plots display the posterior distribution over models (and related metrics).

![Fig1](https://github.com/maheump/Emergence/blob/initialdev/Ideal%20observer/ToyExamples/figs/Emergence_IO_ToyExampleFullIO.jpeg)

The second group of plots display the posterior distribution over models' parameters.

![Fig2](https://github.com/maheump/Emergence/blob/initialdev/Ideal%20observer/ToyExamples/figs/Emergence_IO_ToyExampleFullIO_fig2.jpeg)

The third group of plots display the posterior distribution over change point's position (and related metrics).

![Fig3](https://github.com/maheump/Emergence/blob/initialdev/Ideal%20observer/ToyExamples/figs/Emergence_IO_ToyExampleFullIO_fig3.jpeg)

The fourth group of plots display the expectations regarding the identity of the next observation (and related metrics).

![Fig4](https://github.com/maheump/Emergence/blob/initialdev/Ideal%20observer/ToyExamples/figs/Emergence_IO_ToyExampleFullIO_fig4.jpeg)
