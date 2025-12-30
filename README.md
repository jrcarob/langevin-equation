### The Langevin equation

Inspired by [@hiroloquy](https://github.com/hiroloquy), I replicated his GNU simulation but in R. 

The Langevin equation (introduced by Paul Langevin in 1908 to model **Brownian motion**) is a stochastic differential equation that describes the motion of a particle subject to both deterministic forces (like friction and external potentials) and random fluctuating forces (modeling thermal noise). 
The classic form for a Brownian particle is:

$$m\frac{d^{2}x}{dt^{2}}=-\gamma\frac{dx}{dt}+F(x)+\xi(t)$$

where:
- $m$ is mass,
- $\gamma$ is the friction coefficient,
- $F(x)$ is any systematic force,
- $\xi(t)$ is Gaussian white noise with zero mean and delta-correlated fluctuations.

The code in this repository shows typical trajectory from a Langevin equation simulation, illustrating random Brownian paths (with and without confining potentials).

In **finance**, the Langevin equation models the stochastic behavior of asset prices and other financial variables, capturing both deterministic trends (e.g., drift due to expected returns) and random fluctuations (e.g., market noise or volatility).

Key applications include:

- **Stock price modeling:** The classic Geometric Brownian Motion (GBM) used in the Black-Scholes option pricing model is essentially a Langevin equation in logarithmic form: $dS=\mu Sdt+\sigma SdW$, where $(S)$ is the stock price,
$mu$ is the drift rate, $\sigma$ is volatility, and $(dW)$ is Wiener noise. This enables derivative pricing, risk assessment, and Monte Carlo simulations of price paths.
- **Interest rates and other variables:** Extensions like mean-reverting processes (e.g., Ornstein-Uhlenbeck, a direct Langevin form) model short-term interest rates (Vasicek model) or stochastic volatility (Heston model).
- **Market crashes and extreme events:** Nonlinear Langevin equations incorporate feedback effects (e.g., herding or imitation among traders), creating potential barriers that separate normal fluctuations from crash regimes.
Early models (e.g., Bouchaud and Cont, 1998) use this to explain power-law tails in return distributions and sudden market collapses.
- **Econophysics and data-driven inference:** Researchers estimate Langevin potentials directly from time series data to identify stable price equilibria, metastability, or non-stationary dynamics in stocks, currencies, or portfolios.
This aids risk management, portfolio optimization, and forecasting under uncertainty.

![image](https://github.com/jrcarob/langevin-equation/blob/master/Langevin_eq.gif)
