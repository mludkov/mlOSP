# mlOSP
Regression Monte Carlo for Optimal Stopping


Description: A suite of regression Monte Carlo algorithms utilizing *Machine Learning for Optimal Stopping Problems* (mlOSP). 

Includes both static and
sequential experimental designs. We implement the original Longstaff-Schwartz
and Tsitsiklis-van Roy algorithms, as well as machine learning
approaches that explicitly specify the underlying experimental
designs. The mlOSP template then allows to mix and match the choice
of the regression method, the experimental design, and the
stochastic simulator. Key solver functions are
**osp.prob.design** (original LSM), **osp.fixed.design** (a variety of
space-filling or user-specified designs, generally assumed to be
batched), **osp.seq.design** (sequential designs using a collection of
pre-specified Expected Improvement Criteria), **osp.seq.batch.design** (
sequential design with adaptive batching) and **osp.tvr** (TvR method). 
Also implements the Bouchard-Warin hierarchical adaptive partitioning with
 linear regression (**osp.design.piecewisebw**). The library currently works
with 10+ regression emulators, see documentation. 

The Bermudan_demo *vignette* provides a short illustration with a 2D Bermudan
basket Put. The two demo R files provide the source code for the 
http://arxiv.org/abs/2012.00729 article and the respective benchmarked solvers.

Work partially supported by NSF-1521743 and NSF-1821240.
