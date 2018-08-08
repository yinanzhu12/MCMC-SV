In this repo, I implement three different MCMC ([Markov chain Monte Carlo](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo)) scheme and use them to fit stochastic volatility model.

The stochastic volatility model is a wide family of [models](https://en.wikipedia.org/wiki/Stochastic_volatility). Here we consider the simplest one described in [many places](https://github.com/yinanzhu12/MCMC-SV/blob/master/reference/Bayesian%20Analysis%20of%20Stochastic%20Volatility%20Models.pdf)

All the three schemes are based on Gibbs sampling and differentiate each other in the way they treat a complicated distribution during sampling. The formalisms are described in [my write-up](https://github.com/yinanzhu12/MCMC-SV/blob/master/document/Theory%20and%20Application%20of%20MCMC%20Algorithm%20on%20Stochastic%20Volatility%20Model.pdf) and reference cited there and the code is in code/functions.R

In code/svol.R we apply the model to S&P 500 data in 2007-2010 

