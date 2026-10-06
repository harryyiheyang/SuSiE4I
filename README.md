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

## Factor Variables in Z

An unordered factor (for example sleep or drinking category) can be passed in
`Z` as indicator columns with the baseline level dropped, and `groupint_ind`
lists which columns form each factor, for example
`groupint_ind = list(Drink = c("dr_1", "dr_2", "dr_3"))` (indices or names; no
column in two groups). The interaction columns of one factor with one
main-effect credible set form ONE candidate in the interaction SuSiE: a group
single effect with prior `N(0, V I / d)` on its `d` columns, so `V` is the total
effect variance of the group. Its joint Bayes factor competes with the other
candidates and the levels inside a group do not compete, so an effect spread
over several levels is pooled into one direction. The interacting levels are
described by the per-level `lfsr` (local false sign rate of that single
effect, reported as is; levels are compared with the dropped baseline) in
`interaction_discoveries`, which lists every level of a
selected group with the group's `PIP` and names both sides (`Group1`, `Term1`,
`Group2`, `Term2`, `Pair`). A group enters the refit as one column (its
posterior direction), so `Pvalue` tests the whole group effect. An ordinal
factor can instead enter `Z` as one score column, and haplotypes belong in `X`.
See `example/example_group_int.R`.

Interaction components that SuSiE did not kill (prior variance above zero)
but that did not form a credible set can still enter the refit. Each one takes
its coverage set at `coverage_nonkilled`, purified by dropping members with
absolute correlation below `min_abs_corr` to the lead; if the purified set
still reaches `coverage_nonkilled`, the component enters the refit (built from
the purified set) and is listed in `interaction_discoveries` with
`InCS = FALSE` (the last column) and its `Coverage`. Components that do not
enter the refit are not listed. `coverage_nonkilled` defaults to the smaller of
the interaction CS coverage and 0.8. This helps sparse level-by-level cells,
which often stay below the coverage needed for a credible set.
Main effects still require a credible set. With `groupint_ind`, each `L_int`
component selects a whole group, so the default `L_int = 5` applies.

## Refit Summaries

Selected credible-set summaries are refit jointly in the outcome model. Optional
non-CS residual summaries may be included only as nuisance refit covariates;
they are not reported as credible sets. This improves the linear predictor used
in subsequent iterations while keeping the reported discoveries tied to SuSiE
credible sets.
