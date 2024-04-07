# Towards Markov Chains

## DNA as a Markov chain

Several packages (e.g. [MarkovChainsHammer.jl](https://github.com/sandreza/MarkovChainHammer.jl), [DiscreteMarkovChains.jl](https://github.com/Maelstrom6/DiscreteMarkovChains.jl), etc.) in the Julia ecosystem have been implemented to work with Markov chains with a *state space* of integers, those could be efficient in many ways, but they are clumsy to work with a specialized biological types as in the `BioJulia` ecosystem. Therefore, in the `GeneFinder` package we dedicated some implementations to work with `BioSequence` types so that we can expand the functionality in an efficient way (see complete [API](https://camilogarciabotero.github.io/GeneFinder.jl/dev/api/)).

One important step towards many gene finding algorithms is to represent a DNA sequence as a Markov chain. In this representation a DNA sequence of a reduced alphabet ``\mathscr{A} = \{A,  C,  G ,  T \}`` is draw as a four-vertex graph, where each letter of ``\mathscr{A}`` is a *state (vertex)* and the edges of the graph represent *transitions* from one nucleotide to another in a sequence (e.g. ``A \rightarrow T`` represent a single nucleotide to nucleotide transition). This is also considered more specifically as a Discrete Markov chain (Axelson-Fisk 2015). The complete set of transitions and states of a DNA sequence of alphabet ``\mathscr{A}``.

More formally a Markov chain is a random process where each state is a random variable ``X_{t}`` where ``t \in T`` is a discrete time in a finite sequence ``T`` and the probability to jump from one state into another is *only dependent of the current state*. Therefore a definition of this *Markov property* is given by:

``` math
\begin{align}
P(X_{t} = j |X_{t−1} = i)
\end{align}
```

where ``i, j \in  \mathscr{A}`` . This property led us to generalize a way to calculate the probability of a sequence ``T`` from a process ``(X_{1}...X_{T})`` where each random variable is a nucleotide from ``\mathscr{A}`` so that:

``` math
\begin{align}
P(X_{1} = i_{1},...,X_{T} = i_{T}) = P(X_{1} = i_{1}) \prod_{t=2}^{T} P(X_{t} = i_{t} | X_{t−1} = i_{t−1})
\end{align}
```

Note that previous equations has two terms, a initial probability ``P(X_{1} = i_{1})`` and the the product of all transitions beginning at ``t = 2``. So, to calculate the initial probability distribution of each of the nucleotides of a string ``T`` with the alphabet 𝒜 we can first calculate the transition probability matrix ``\widehat{\mathscr{M}}`` out of the frequency count of the transitions. In an alphabet 𝒜 we got ``4^2`` transitions of one order, that is the ``AA, AC, AG, ...`` which coincides with the frequency of the dinucleotides in the sequence. So we can later in fact build a ``4 x 4`` matrix representing all the transitions. For instance in a DNA sequence ``T`` of 24 nucleotides:  

``` math
CCTCCCGGACCCTGGGCTCGGGAC
```

We can calculate each frequency nucleotide to any other nucleotide ``\widehat{m}_{ij} = \frac{c_{ij}}{c_{i}}`` where ``c_{ij}`` is the actual count of the dinucleotide, and therefore ``c_{i}`` is the counts of the nucleotide ``i`` to any other nucleotide and build the transition probability matrix:

``` math
\begin{bmatrix}
   & \text{A} & \text{C} & \text{G} & \text{T} \\
\text{A} & 0.00 & 1.00 & 0.00 & 0.00 \\
\text{C} & 0.00 & 0.56 & 0.22 & 0.30 \\
\text{G} & 0.25 & 0.12 & 0.62 & 0.00 \\
\text{T} & 0.00 & 0.67 & 0.33 & 0.00 \\
\end{bmatrix}
```

It is noteworthy that initial probabilities can also be obtained from the counts of each nucleotide transitions ``c_{ij}`` over the total sum of the dinucleotide counts ``c_{k}``:

``` math
\begin{align}
\widehat{\pi}_{i} = \frac{c_{i}}{\sum_{k}c_{k}}
\end{align}
```

That way for the previous example example we can can calculate the initial probabilities ``\widehat{\pi} = (0.08,0.43,0.34,0.13)``. Both set of probabilities composed a *transition model* that can be used to predict the probability of any DNA sequence using equation (2).

## References

Axelson-Fisk, Marina. 2015. *Comparative Gene Finding*. Vol. 20. Computational Biology. London: Springer London. <http://link.springer.com/10.1007/978-1-4471-6693-1>.