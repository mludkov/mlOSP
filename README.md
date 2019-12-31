# mlOSP
Regression Monte Carlo for Optimal Stopping


Description: A suite of regression Monte Carlo algorithms utilizing *Machine Learning for Optimal Stopping Problems* (mlOSP). 

Includes both static and
sequential experimental designs. We implement the original Longstaff-Schwartz
and Tsitsiklis-van Roy algorithms, as well as machine learning
approaches that explicitly capture the underlying experimental
designs. The mlOSP template then allows to mix and match the choice
of the regression method, the experimental design and the
stochastic simulator. The library directly accepts function hooks
for the option payoff and the path generation. Key functions are
**osp.prob.design** (original LSM), **osp.fixed.design** (a variety of
space-filling or user-specified designs, generally assumed to be
batched), **osp.seq.design** (sequential designs using a collection of
pre-specified Expected Improvement Criteria). Also implements the
Bouchard-Warin adaptive partitioning with linear regression
(**osp.design.piecewisebw**). 

Work partially supported by NSF-1521743.
