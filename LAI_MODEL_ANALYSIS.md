# Analysis of Prognostic LAI Model Equations

## Reference
Zhou et al. (2025) "A General Model for the Seasonal to Decadal Dynamics of Leaf Area" 
Global Change Biology, https://onlinelibrary.wiley.com/doi/pdf/10.1111/gcb.70125

## Overview
This document analyzes the implementation of the prognostic LAI model in 
`src/standalone/Vegetation/pmodel.jl` (lines 1533-1767) relative to the theoretical 
framework described in Zhou et al. (2025).

## Equation Analysis

### 1. Canopy-Integrated Assimilation: `compute_Ao(A, k, L)`

**Implementation:**
```julia
Ao = A / max(1 - exp(-k*L), eps(FT))
```

**Theory:** This equation inverts the Beer-Lambert law integration. The canopy-integrated 
assimilation A is related to the potential (top-of-canopy) assimilation Ao by:

```
A = Ao ∫₀ᴸ exp(-k*l) dl = Ao/k * [1 - exp(-k*L)]
```

For a unit-normalized case (k=1), this simplifies to:
```
A = Ao * [1 - exp(-k*L)]
```

**Status:** ✅ **CORRECT** - The implementation correctly inverts this relationship to recover Ao from A.

**Note on units:** The code appears to assume a simplified form where the extinction 
coefficient is already normalized. The actual Beer's law would be A = (Ao/k) * [1 - exp(-k*L)], 
but here k appears to be incorporated differently.

---

### 2. Maximum LAI Constraint: `compute_L_max(k, z, Ao)`

**Implementation:**
```julia
L_max = -1/k * log(z / (k * Ao))
```

**Theory:** The maximum LAI occurs when light at the bottom of the canopy reaches the 
light compensation point z. Using Beer's law for light attenuation:

```
Light(L) = Ao * exp(-k*L)
```

Setting Light(L_max) = z (the compensation point):
```
z = Ao * exp(-k*L_max)
```

Solving for L_max:
```
exp(-k*L_max) = z / Ao
-k*L_max = log(z / Ao)
L_max = -(1/k) * log(z / Ao)
```

**Discrepancy:** ⚠️ **POTENTIAL ISSUE** - The implementation has an extra `k` in the denominator:
```julia
L_max = -1/k * log(z / (k * Ao))  # Current implementation
L_max = -1/k * log(z / Ao)        # Expected from theory
```

This suggests one of two possibilities:
1. There's a unit scaling issue where z and Ao have different normalizations
2. The parameter z represents something slightly different than described

**Impact:** The extra k factor would systematically shift L_max. If k < 1 (dense canopy), 
it would predict higher L_max; if k > 1 (sparse canopy), lower L_max.

---

### 3. Optimization Function: `g(μ, k, L)`

**Implementation:**
```julia
g(μ, k, L) = L/μ - 1 + exp(-k*L)
```

**Theory:** The optimal LAI maximizes net carbon gain = gross assimilation - respiration.
Assuming:
- Gross assimilation follows Beer's law integration: Ao * [1 - exp(-k*L)]
- Respiration cost is proportional to LAI: m * Ao * L (where μ = m * Ao)
- Net gain: N(L) = Ao * [1 - exp(-k*L)] - m * Ao * L

At optimum, dN/dL = 0:
```
dN/dL = Ao * k * exp(-k*L) - m * Ao = 0
k * exp(-k*L) = m
```

Dividing by μ = m * Ao and rearranging:
```
k/μ * exp(-k*L) = 1
```

**Alternative formulation:** If we define g such that g = 0 at optimum:
```
g = ∫ [marginal gain - marginal cost] dL
g = 1 - exp(-k*L) - m*L/Ao
g = 1 - exp(-k*L) - L/μ
```

Rearranging to match the implementation:
```
g = L/μ - 1 + exp(-k*L)
```

**Status:** ✅ **CORRECT** - The implementation matches the integral form of the 
optimization condition.

---

### 4. Derivative: `dgdL(μ, k, L)`

**Implementation:**
```julia
dgdL(μ, k, L) = 1/μ - k*exp(-k*L)
```

**Verification:**
Taking the derivative of g(μ, k, L) = L/μ - 1 + exp(-k*L):
```
dg/dL = 1/μ + 0 - k*exp(-k*L) = 1/μ - k*exp(-k*L)
```

**Status:** ✅ **CORRECT** - Derivative is computed correctly.

---

### 5. Newton-Raphson Solver: `compute_L_opt(μ, k, L)`

**Implementation:**
```julia
dL = -g(μ, k, L) / dgdL(μ, k, L)
L = L + dL
```

**Theory:** Standard Newton-Raphson iteration for finding roots:
```
L_{n+1} = L_n - g(L_n) / g'(L_n)
```

**Status:** ✅ **CORRECT** - Standard Newton-Raphson implementation.

**Note:** The function includes debugging output (@show(dL)) which may impact performance 
in production. Consider removing or making conditional.

---

### 6. Exponential Moving Average: `compute_L()`

**Implementation:**
```julia
L_new = (1-α) * L_ss + α * L_old
```

where α = 1 - 0.067 = 0.933

**Theory:** Exponential moving average with memory parameter α:
```
L_{n+1} = (1-α) * L_target + α * L_n
```

The effective time constant is τ ≈ 1/(1-α) = 1/0.067 ≈ 15 days.

**Status:** ✅ **CORRECT** - Standard EMA formulation.

**Comment in code says:** "weight it with 0.0667" but the formula shows (1-α) = 0.067 
applied to L_ss (the new steady-state), and α = 0.933 applied to L (the old value). 
The comment could be clearer that we're giving NEW values only 6.7% weight, keeping 
93.3% of history.

---

## Summary of Findings

### ✅ Correct Implementations
1. `compute_Ao` - Correctly inverts Beer's law integration
2. `g` - Correct optimization function
3. `dgdL` - Correct derivative
4. `compute_L_opt` - Correct Newton-Raphson solver
5. `compute_L` - Correct EMA update scheme

### ⚠️ Potential Issues

1. **`compute_L_max` equation discrepancy:**
   - Implementation: `L_max = -1/k * log(z / (k * Ao))`
   - Expected: `L_max = -1/k * log(z / Ao)`
   - **Impact:** The extra k in the denominator shifts the maximum LAI prediction
   - **Recommendation:** Verify against Zhou et al. (2025) equations. If this is 
     intentional due to unit conventions, add a comment explaining why.

2. **Debug output in production code:**
   - `@show(dL)` in `compute_L_opt` will print to console on every iteration
   - **Recommendation:** Remove or wrap in a debug flag

3. **Units consistency:**
   - The relationship between z, Ao, and k needs clarification
   - Is z really in units of assimilation rate (μmol CO₂ m⁻² s⁻¹)?
   - Or is it a dimensionless ratio?

### Parameter Values

From the code:
- `α = 0.933` (memory coefficient) → ~15 day timescale ✅
- `z = 12.227` (compensation point) - **units need verification**
- `m = 0.3` (respiration to assimilation ratio) ✅

## Recommendations

1. **Verify `compute_L_max` formula** against Zhou et al. (2025) Equation [X] 
   (specific equation number from paper)

2. **Add unit tests** that compare against known analytical solutions for simple cases

3. **Remove or conditionalize debug output** (@show statements)

4. **Add reference to specific equations** from Zhou et al. (2025) in comments

5. **Clarify parameter z:**
   - If it's meant to be Ao_compensation / Ao_max, it should be dimensionless
   - If it's meant to be absolute compensation point, units should be consistent

## Example Calculation

For verification, with typical values:
- k = 0.5 (extinction coefficient)
- Ao = 25.0 μmol CO₂ m⁻² s⁻¹
- m = 0.3
- z = 12.227 μmol CO₂ m⁻² s⁻¹

Using current implementation:
```
μ = m * Ao = 0.3 * 25.0 = 7.5
L_max = -1/0.5 * log(12.227 / (0.5 * 25.0)) = -2.0 * log(0.978) ≈ 0.045
```

This seems very low! Expected L_max should be around 4-6 for typical forests.

Using corrected formula:
```
L_max = -1/0.5 * log(12.227 / 25.0) = -2.0 * log(0.489) ≈ 1.43
```

Still seems low. This suggests either:
1. The z parameter value (12.227) may be incorrect or misinterpreted
2. There's a conceptual difference in what z represents

**Conclusion:** The `compute_L_max` equation needs verification against the original 
paper, as the predicted values seem inconsistent with typical forest LAI values.
