% PNLSS 1.0
% A polynomial nonlinear state-space toolbox for Matlab
% September 2, 2016
%
% TUTORIAL
%   sTutorialSISO                       - Tutorial that shows how to simulate and estimate a SISO PNLSS model.
%
% SIGNAL GENERATION AND TRANSIENT HANDLING
%   fComputeIndicesTransient           - Computes indices for transient handling for periodic signals before filtering.
%   fComputeIndicesTransientRemoval    - Computes indices for transient handling for periodic signals after filtering.
%   fComputeIndicesTransientRemovalArb - Remove transients from arbitrary data.
%   fMultisine                         - Generate a random-phase multisine signal with a flat spectrum.
%
% NONPARAMETRIC FREQUENCY RESPONSE ESTIMATION
%   fCovarFrf                          - Computes frequency response matrix and noise and total covariance matrix from input/output spectra.
%   fCovarY                            - Compute covariance matrix output spectra due to output noise from time domain output signals.
%
% PARAMETRIC SUBSPACE IDENTIFICATION
%   fFreqDomSubSpace                   - Estimate state-space model from Frequency Response Function (or Matrix).
%   fFreqDomSubSpaceCT                 - Estimate state-space model from Frequency Response Function (or Matrix) (continuous-time).
%   fIsUnstable                        - Determines if a linear state-space model is unstable.
%   fJacobFreqSS                       - Compute Jacobians of the unweighted errors w.r.t. elements A, B, and C matrices.
%   fLevMarqFreqSSz                    - Optimize state-space matrices using Levenberg-Marquardt.
%   fLoopSubSpace                      - Loop frequency-domain subspace method over multiple model orders and sizes of the extended observability matrix.
%   fss2frf                            - Compute frequency response function from state-space parameters (discrete-time).
%   fss2frfCT                          - Compute frequency response function from state-space parameters (continuous-time).
%   fStabilize                         - Stabilize a linear state-space model.
%   fWeightJacobSubSpace               - Adds weighting to an unweighted Jacobian.
%
% NONLINEAR OPTIMIZATION
%   fLMnlssWeighted                    - Optimize PNLSS model using weighted Levenberg-Marquardt algorithm.
%   fLMnlssWeighted_x0u0               - Optimize PNLSS model and initial conditions using weighted Levenberg-Marquardt algorithm.
%   fComputeJF                         - Compute Jacobian of F*eta w.r.t. active terms in eta.
%   fEdwdx                             - Multiply a matrix E with the derivative w.r.t. x of a polynomial w(x,u).
%   fEdwdu                             - Multiply a matrix E with the derivative w.r.t. u of a polynomial w(x,u).
%   fJNL                               - Compute Jacobian w.r.t. A, B, and E by filtering an alternative state-space model.
%   fJx0                               - Compute Jacobian w.r.t. x0 by filtering an alternative state-space model.
%   fJu0                               - Compute Jacobian w.r.t. u0 by filtering an alternative state-space model.
%   sJacobianAnalytical                - Script to compute the Jacobians in a nonlinear state-space model.
%   sJacobianAnalytical_x0u0           - Script to compute the Jacobians in a nonlinear state-space model.
%
% MODEL CONSTRUCTION AND SIMULATION
%   fCreateNLSSmodel                   - Create polynomial nonlinear state-space model from initial linear state-space model.
%   fSelectActive                      - Select active elements in E or F matrix.
%   fSScheckDims                       - Check consistency of state-space dimensions.
%   fFilterNLSS                        - Calculate the output and the states of a nonlinear state-space model with transient handling.
%   fFilterspeedNL                     - Calculate the output and the states of a nonlinear state-space model without transient handling.
%
% UTILITY
%   fHerm                              - Average of square matrix and its Hermitian transpose.
%   fMetricPrefix                      - Returns an appropriate metric prefix and a corresponding scaling factor.
%   fNormalizeColumns                  - Normalizes the columns of a matrix with their rms value.
%   fOne                               - Constructs a matrix with only one one, and zeros elsewhere.
%   fPlotFrfMIMO                       - Make amplitude versus frequency plots of the elements of a frequency response matrix.
%   fReIm                              - Stacks the real and imaginary part of a matrix on top of each other.
%   fSqrtInverse                       - Computes B, such that B*B = inv(A).
%   fVec                               - Vectorization of a matrix or tensor.
%   fCombinations                      - Lists all nonlinear terms in a multivariate polynomial.
%   fTermNL                            - Construct polynomial terms.
%
%
%	Copyright (c) Vrije Universiteit Brussel – dept. ELEC
%   All rights reserved.
%   Software can be used freely for non-commercial applications only.
%   Disclaimer: This software is provided “as is” without any warranty.
