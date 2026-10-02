# SuSiE4I

SuSiE4I implements iterative SuSiE fine-mapping for main effects and
interaction effects under Gaussian, generalized linear, binary, and Cox
outcomes. The implementation alternates between constructing local quadratic
working problems, running `susieR::susie_ss()` on blockwise sufficient
statistics, and refitting the selected credible-set summaries in the outcome
model.

## GLM and Extended GLM Path

For a GLM or an extended GLM family, the mean and variance models are

```math
g(\mu_i) = \eta_i, \qquad
\mathrm{Var}(Y_i \mid \mu_i) = \phi V(\mu_i).
```

At the current fit, the code obtains the usual IRLS working response and
weights,

```math
z_i =
\eta_i +
\frac{y_i - \mu_i}{\partial \mu_i / \partial \eta_i},
\qquad
\omega_i =
\frac{(\partial \mu_i / \partial \eta_i)^2}{V(\mu_i)}.
```

After whitening by the square root of the IRLS weights, the local working model
has the form

```math
z_i^\ast =
\eta_0^\ast + X_i^\ast \beta + W_i^\ast \gamma + \varepsilon_i,
\qquad
\varepsilon_i \sim N(0,\sigma^2).
```

For the canonical IRLS quadratic approximation, the theoretical working
residual variance is 1. In finite samples, however, the first few
outer iterations may use imperfect estimates of the linear predictor and, for
extended families, imperfect estimates of dispersion or family parameters.
Therefore SuSiE4I allows `susie_ss()` to estimate

```math
\sigma^2 \in [0.1, 1.01],
```

initialized at `residual_variance = 0.5`. This keeps the working likelihood
close to the unit-variance IRLS target while preventing occasional overly large
variance estimates from unnecessarily reducing power.

The GLM branch is used for families supported by the `mgcv` IRLS machinery,
including Poisson, negative binomial, Tweedie, beta regression, and scaled
t-type responses. At each outer iteration SuSiE4I:

1. fits the current `mgcv::gam()` or `mgcv::bam()` model;
2. extracts the working response and weights;
3. builds blockwise sufficient statistics for main effects;
4. constructs the selected interaction/environment design;
5. builds blockwise sufficient statistics for interaction effects;
6. refits the selected credible-set summaries and updates the linear predictor.

Large matrix products are evaluated through blockwise routines so that the
method does not need to explicitly form all dense intermediate projection
matrices.

## Binary Outcomes

Binary outcomes use the standard GLM IRLS construction above with
`family = binomial(link = "logit")`. The working response, weights, sufficient
statistics, and final refit all remain on the conventional logistic IRLS scale.

## Cox Outcomes

For right-censored survival outcomes, SuSiE4I uses a Cox score path based on
the partial likelihood. The Cox model does not provide an observed Gaussian
response in the same way as the GLM working-response construction. Instead,
the code constructs score and observed-information sufficient statistics from
the risk sets,

```math
X^\top I X, \qquad X^\top U,
```

where `U` is the Cox score contribution and `I` is the observed information
under the current linear predictor.

Because there is no explicit working response with a fixed unit residual
variance, the Cox path uses

```math
y^\top y = n - 1
```

as an information-scale normalization and leaves a small degree of freedom for
`susie_ss()` to estimate the residual variance. In contrast to the GLM and
extended GLM IRLS paths, this variance is not theoretically forced to equal
one. The default Cox setting is therefore

```math
\sigma^2 \in [0.1, 1.01],
```

with `residual_variance = 0.5` and
`estimate_residual_variance = TRUE`.

The Cox branch uses the same high-level iteration as the other branches:
score-based sufficient statistics, SuSiE main-effect fitting, interaction
construction from selected credible sets, SuSiE interaction fitting, and a
final Cox partial-likelihood refit on the selected summaries.

## Grouped Factors in Z (Haplotypes)

Factor-coded covariates such as haplotypes can be passed in `Z` as reference
coded indicator columns together with `z_groups`, for example
`z_groups = list(HapA = c("HA_a1", "HA_a2"), HapB = c("HB_b1"))`. The factors
stay unpenalized in `Z`. Instead of treating each indicator as an unrelated
covariate, the interaction generator then forms one column per level times
each main credible-set column and one column per level-by-level product of two
different factors. Products within a factor are never formed. A column is
dropped when fewer than `min_group_int_obs` observations (default 100) carry
the level, or both levels for a level-by-level product. Selection in the
interaction stage is unchanged (`susieR::susie_ss` on these columns), and
`interaction_discoveries` names the factor and level on each side of a
selected interaction (`Factor1`, `Level1`, `Factor2`, `Level2`, `Pair`), so
the specific interacting haplotype levels can be read off directly. See
`example/example_group_int.R`.

## Refit Summaries

Selected credible-set summaries are refit jointly in the outcome model. Optional
non-CS residual summaries may be included only as nuisance refit covariates;
they are not reported as credible sets. This improves the linear predictor used
in subsequent iterations while keeping the reported discoveries tied to SuSiE
credible sets.
