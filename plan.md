# Form Factor Refactoring Plan

## Overview
Refactor Blatt-Weisskopf form factors from lineshapes into VertexFunction using HadronicLineshapes.jl, maintaining exact numerical compatibility through incremental changes.

## Current State
- Original lineshape: `1 / (m^2 - σ - 1im * m * Γ) * (p / p0)^l * (q / q0)^minL * sqrt(F²(l, p, p0, dR) * F²(minL, q, q0, dΛc))`
- Form factors are embedded in the lineshape calculation
- Tests pass with exact numerical values

## Target State
- Lineshape: `1 / (m^2 - σ - 1im * m * Γ) * (q / q0)^minL * sqrt(F²(minL, q, q0, dΛc))` (keep q-dependent terms)
- Hij vertex: `BlattWeisskopf{l}(dR)` form factor
- HRk vertex: `NoFormFactor()` (q-dependent terms stay in lineshape)
- Mathematical equivalence: `HRk * X_new * Hij = HRk_old * X_old * Hij_old`

## Mathematical Decomposition

### Original Expression
```
1 / (m^2 - σ - 1im * m * Γ) * (p / p0)^l * (q / q0)^minL * sqrt(F²(l, p, p0, dR) * F²(minL, q, q0, dΛc))
```

### Target Decomposition
```
FF_Rk * X * FF_ij
```

Where:
- `FF_Rk = BlattWeisskopf{minL}(dΛc)` (HRk vertex form factor)
- `FF_ij = BlattWeisskopf{l}(dR)` (Hij vertex form factor)
- `X = 1 / (m^2 - σ - 1im * m * Γ) * (1 / dR / p0)^l * (1 / dΛc / q0)^minL * sqrt(F²(l, 0, p0, dR) * F²(minL, 0, q0, dΛc) * factor(l) * factor(minL))`

With:
- `factor(l) = l == 2 ? 1/9 : 1`
- `factor(minL) = minL == 2 ? 1/9 : 1`

### Validation Requirements
**CRITICAL**: This mathematical decomposition must be validated separately and upfront before any code modifications:

1. **Mathematical equivalence proof**: Verify that `FF_Rk * X * FF_ij = Original_Expression`
2. **Numerical validation**: Test with multiple kinematic points
3. **Edge case testing**: Test with l=0, l=1, l=2, minL=0, minL=1, minL=2
4. **Form factor interface validation**: Ensure HadronicLineshapes.jl `BlattWeisskopf{L}(d)` matches expected behavior

### Pre-Implementation Validation Steps
1. **Create validation script** (`validate_decomposition.jl`)
2. **Test mathematical equivalence** with symbolic or high-precision arithmetic
3. **Verify form factor behavior**:
   - `BlattWeisskopf{l}(dR)(p, p0, σ)` should give `(p * dR)^l`
   - `BlattWeisskopf{minL}(dΛc)(q, q0, σ)` should give `(q * dΛc)^minL`
4. **Test compensation factors**:
   - Verify `(1 / dR / p0)^l * (1 / dΛc / q0)^minL` correctly compensates
   - Verify `sqrt(F²(l, 0, p0, dR) * F²(minL, 0, q0, dΛc))` handles reference values
   - Verify `factor(l) * factor(minL)` handles the l=2 case correctly

## Step-by-Step Implementation Plan

### Phase 1: Setup and Preparation
**Goal**: Establish baseline and prepare infrastructure

1. **Create backup branch**
   ```bash
   git checkout -b form-factor-refactoring
   git add -A && git commit -m "Baseline before form factor refactoring"
   ```

2. **Add HadronicLineshapes.jl dependency**
   - Add to `Project.toml`
   - Add `using HadronicLineshapes` to main module
   - Test that package loads correctly

3. **Create test harness**
   - Create `test_form_factors.jl` with exact numerical checks
   - Store reference values for key test cases
   - Ensure all tests pass before starting

### Phase 2: Mathematical Validation
**Goal**: Validate the mathematical decomposition before any code changes

4. **Create validation script** (`validate_decomposition.jl`)
   - Implement original expression
   - Implement target decomposition
   - Test mathematical equivalence with high precision

5. **Validate form factor behavior**
   - Test `BlattWeisskopf{L}(d)` with known inputs
   - Verify `BlattWeisskopf{l}(dR)(p, p0, σ) = (p * dR)^l`
   - Verify `BlattWeisskopf{minL}(dΛc)(q, q0, σ) = (q * dΛc)^minL`

6. **Test compensation factors**
   - Verify `(1 / dR / p0)^l * (1 / dΛc / q0)^minL` compensation
   - Test `F²(l, 0, p0, dR)` and `F²(minL, 0, q0, dΛc)` at reference values
   - Test `factor(l) * factor(minL)` for l=2, minL=2 cases

7. **Edge case validation**
   - Test with l=0, l=1, l=2
   - Test with minL=0, minL=1, minL=2
   - Test with various kinematic points
   - Ensure numerical stability

### Phase 3: Interface Compatibility
**Goal**: Understand how to use HadronicLineshapes.jl correctly

8. **Study HadronicLineshapes.jl API**
   - Test `BlattWeisskopf{L}(d)` constructor
   - Understand calling convention: `ff(m0sq, m1sq, m2sq)` vs `ff(p, p0, σ)`
   - Create adapter if needed

9. **Create compatibility layer**
   - Implement `BlattWeisskopfAdapter{L}` that bridges ThreeBodyDecays interface
   - Test with simple cases to ensure mathematical equivalence

### Phase 4: Incremental Form Factor Migration
**Goal**: Move form factors one by one, maintaining compatibility

10. **Step 1: Move Hij form factor (p-dependent terms)**
   - **Before**: Lineshape has `(p / p0)^l * sqrt(F²(l, p, p0, dR))`
   - **After**: 
     - Lineshape: Remove `(p / p0)^l * sqrt(F²(l, p, p0, dR))`
     - Hij vertex: Add `BlattWeisskopf{l}(dR)`
     - Add compensation: `1 / ff_hij(p0, p0, σ)` in lineshape
   - **Test**: Verify exact numerical compatibility

11. **Step 2: Move HRk form factor (q-dependent terms)**
   - **Before**: Lineshape has `(q / q0)^minL * sqrt(F²(minL, q, q0, dΛc))`
   - **After**:
     - Lineshape: Remove `(q / q0)^minL * sqrt(F²(minL, q, q0, dΛc))`
     - HRk vertex: Add `BlattWeisskopf{minL}(dΛc)`
     - Add compensation: `1 / ff_hrk(q0, q0, σ)` in lineshape
   - **Test**: Verify exact numerical compatibility

### Phase 5: Cleanup and Optimization
**Goal**: Remove temporary code and optimize

12. **Remove compensation factors**
   - Once form factors are properly distributed, remove compensation
   - Verify that `HRk * X * Hij` gives same result as original

13. **Clean up code**
   - Remove `BlattWeisskopfWrapper` if no longer needed
   - Update documentation
   - Run full test suite

## Critical Warnings and Considerations

### Mathematical Decomposition Validation
- **CRITICAL**: The mathematical decomposition must be validated BEFORE any code changes
- The expression `FF_Rk * X * FF_ij = Original_Expression` must be proven mathematically
- The compensation factors `(1 / dR / p0)^l * (1 / dΛc / q0)^minL` are non-trivial and must be tested
- The `factor(l) * factor(minL)` terms for l=2, minL=2 cases are critical and must be validated
- **DO NOT PROCEED** without passing the mathematical validation phase

### Mathematical Equivalence
- **CRITICAL**: `HRk_new * X_new * Hij_new = HRk_old * X_old * Hij_old`
- Form factors must be applied at the correct kinematic points
- Compensation factors must be calculated at reference momenta

### Interface Compatibility
- **WARNING**: HadronicLineshapes.jl `BlattWeisskopf{L}(d)` expects `(m0sq, m1sq, m2sq)`
- ThreeBodyDecays expects `(p, p0, σ)` interface
- May need adapter layer

### Numerical Stability
- **WARNING**: Form factors can cause numerical instabilities
- Test with edge cases (low momentum, threshold regions)
- Ensure proper handling of complex arguments

### Testing Strategy
- **CRITICAL**: Test after each step
- Store reference values before starting
- Use exact numerical comparisons, not approximate
- Test multiple kinematic points

## Rollback Plan
If any step fails:
1. Revert to previous working state
2. Analyze the failure
3. Adjust approach
4. Try again with smaller steps

## Success Criteria
- [ ] All tests pass with exact numerical compatibility
- [ ] Form factors are properly distributed between lineshape and vertices
- [ ] Code is cleaner and more modular
- [ ] Documentation is updated
- [ ] Performance is maintained or improved

## File Modifications Required

### Core Files
- `src/lineshapes.jl`: Remove form factors, add compensation
- `src/io.jl`: Update Hij vertex creation
- `src/mapping.jl`: Update HRk vertex creation
- `src/Lc2ppiKSemileptonicModelLHCb.jl`: Add HadronicLineshapes import

### New Files
- `test_form_factors.jl`: Comprehensive form factor testing
- `src/form_factor_adapter.jl`: Adapter for HadronicLineshapes.jl (if needed)

### Configuration
- `Project.toml`: Add HadronicLineshapes dependency
- `Manifest.toml`: Will be updated automatically

## Timeline Estimate
- Phase 1: 1-2 hours
- Phase 2: 2-3 hours  
- Phase 3: 4-6 hours
- Phase 4: 1-2 hours
- **Total**: 8-13 hours

## Risk Mitigation
- Frequent commits at each milestone
- Comprehensive testing at each step
- Rollback capability
- Incremental approach reduces risk of major failures
