---
# Fill in the fields below to create a basic custom agent for your repository.
# The Copilot CLI can be used for local testing: https://gh.io/customagents/cli
# To make this agent available, merge this file into the default repository branch.
# For format details, see: https://gh.io/customagents/config

name: grid/cell migration
description:
---

# AI Agent Migration Guide: Cell and Grid Modernization

## Mission Overview

You are migrating a C++ coastal modeling codebase from old-style class-based cell access to modern struct-based direct access. Your task is to:

1. **Replace old grid/cell access patterns** with modern equivalents
2. **Verify attributes exist** in the new Cell structure
3. **Add simple missing attributes** to Cell struct when straightforward
4. **Flag complex cases** with `// TODO AI-MIGRATION:` comments for manual review
5. **Maintain code functionality** - no behavior changes, only syntax updates

## Context

### Old System (Being Replaced)
- **Grid access:** `m_pRasterGrid->m_Cell[nX][nY]`
- **Cell methods:** `pCell->dGetWaveHeight()`, `pCell->SetWaveHeight(5.0)`
- **Component access:** `pCell->pGetLandform()->nGetLFCategory()`
- **Old classes:** `CGeomCell`, `CRWCellLandform`, `CRWCellLayer`, `CRWCellSediment`, `CRWCellTalus`

### New System (Target)
- **Grid access:** `m_pRasterGrid->Grid()(nX, nY)`
- **Cell struct:** Direct field access `cell.wave_height`, `cell.wave_angle`
- **Nested components:** `cell.landform.category`, `cell.layers[0].unconsolidated.sand`
- **Main file:** `src/cell_struct.h` contains complete Cell definition

---

## Critical Files

### Read These First

1. **`src/cell_struct.h`** - Complete Cell structure definition
   - Contains all field names (NO Hungarian notation)
   - Has nested structs: `Talus`, `Sediment`, `Layer`, `Landform`
   - Review this to know what fields exist

2. **`src/member_mapping.txt`** - Old → New name mappings
   - Format: `m_dWaveHeight -> wave_height`
   - Use this for name translation

3. **`src/MIGRATION_CHEATSHEET.md`** - Quick reference patterns

---

## Decision Tree

For each old access pattern you encounter:

```
1. Is it a simple field access (getter/setter)?
   YES → Apply Pattern 1 (Simple Field Access)
   NO → Continue

2. Is it accessing components (landform/layer/sediment/talus)?
   YES → Apply Pattern 2 (Nested Component Access)
   NO → Continue

3. Is it a method call (calculation/logic)?
   YES → Check if method exists in Cell struct
       - If method is simple → Apply Pattern 3 (Method Migration)
       - If method is complex → FLAG with TODO comment
   NO → Continue

4. Is the field missing from Cell struct?
   YES → Check if it's simple (bool/int/double with clear purpose)
       - If simple → Apply Pattern 4 (Add Missing Field)
       - If complex → FLAG with TODO comment
   NO → FLAG with TODO comment (unknown case)
```

---

## Migration Patterns

### Pattern 1: Simple Field Access (90% of cases)

#### Input Recognition
```cpp
// Old getter patterns:
pCell->dGetWaveHeight()
pCell->dGetWaveAngle()
pCell->dGetSeaDepth()
pCell->bIsInContiguousSea()
pCell->nGetCoastlineID()

// Old setter patterns:
pCell->SetWaveHeight(5.0)
pCell->SetInContiguousSea(true)
pCell->SetCoastlineID(42)
```

#### Name Translation Rules
1. Remove `dGet`, `bIs`, `nGet`, `Set` prefixes
2. Convert CamelCase to snake_case
3. Check `member_mapping.txt` for confirmation

**Examples:**
- `dGetWaveHeight()` → `wave_height`
- `bIsInContiguousSea()` → `in_contiguous_sea`
- `nGetCoastlineID()` → `coastline_id`
- `SetWaveHeight(val)` → `wave_height = val`

#### Output Pattern
```cpp
// BEFORE:
m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight()

// AFTER:
m_pRasterGrid->Grid()(nX, nY).wave_height

// Or better (cache reference):
Cell& cell = m_pRasterGrid->Grid()(nX, nY);
cell.wave_height
```

#### Multi-Access Optimization
```cpp
// BEFORE (multiple accesses):
double h = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();
double a = m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle();
m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth(10.0);

// AFTER (cache cell reference):
Cell& cell = m_pRasterGrid->Grid()(nX, nY);
double h = cell.wave_height;
double a = cell.wave_angle;
cell.sea_depth = 10.0;
```

### Pattern 2: Nested Component Access

#### 2A: Landform Access
```cpp
// BEFORE:
pCell->pGetLandform()->nGetLFCategory()
pCell->pGetLandform()->SetAccumWaveEnergy(100.0)
pCell->pGetLandform()->dGetCliffNotchApexElev()

// AFTER:
cell.landform.category
cell.landform.accum_wave_energy = 100.0
cell.landform.cliff.notch_apex_elev
```

**Landform fields in Cell struct:**
```cpp
struct Landform {
   int category;
   int coast;
   int point_on_coastline;
   double accum_wave_energy;

   struct Cliff {
      double notch_apex_elev;
      double notch_incision;
   } cliff;
};
```

#### 2B: Layer Access
```cpp
// BEFORE:
pCell->pGetLayerAboveBasement(n)->dGetTotalThickness()
pCell->pGetLayerAboveBasement(0)->pGetUnconsolidatedSediment()->dGetSandDepth()

// AFTER:
cell.layers[n].total_thickness()
cell.layers[0].unconsolidated.sand
```

**Layer structure:**
```cpp
struct Layer {
   Sediment unconsolidated;
   Sediment consolidated;
   std::unique_ptr<Talus> talus;

   bool has_talus() const noexcept;
   double total_thickness() const noexcept;
   double uncons_depth() const noexcept;
   double cons_depth() const noexcept;
};
```

#### 2C: Sediment Access
```cpp
// BEFORE:
pLayer->pGetUnconsolidatedSediment()->dGetSandDepth()
pLayer->pGetConsolidatedSediment()->SetFineDepth(5.0)

// AFTER:
layer.unconsolidated.sand
layer.consolidated.fine = 5.0
```

**Sediment fields:**
```cpp
struct Sediment {
   double fine;
   double sand;
   double coarse;
   double notch_fine_lost;
   double notch_sand_lost;
   double notch_coarse_lost;
   // ... more fields

   double total_depth() const noexcept;  // Method
};
```

#### 2D: Talus Access (Optional Member)
```cpp
// BEFORE:
if (pLayer->bHasTalus())
{
   pTalus = pLayer->pGetTalus();
   double sand = pTalus->dGetSandDepth();
}

// AFTER:
if (layer.has_talus())
{
   double sand = layer.talus->sand;
}

// Or create if needed:
auto& talus = layer.get_or_create_talus();
talus.sand = 1.0;
```

### Pattern 3: Method Migration

#### 3A: Methods That Exist in Cell Struct
These methods are preserved as inline helpers:

```cpp
// Layer methods (keep as-is):
layer.total_thickness()
layer.has_talus()
layer.uncons_depth()
layer.cons_depth()

// Sediment methods:
sediment.total_depth()
```

#### 3B: Methods to Convert to Direct Calculation
```cpp
// BEFORE:
pCell->nGetNumLayers()

// AFTER:
cell.layers.size()

// BEFORE:
pCell->pGetLayerAboveBasement(n)

// AFTER:
cell.layers[n]
```

### Pattern 4: Add Missing Field

#### When to Add
Add a field if ALL conditions are met:
1. ✅ Field name is clear and descriptive
2. ✅ Type is simple: `bool`, `int`, or `double`
3. ✅ Purpose is obvious from context
4. ✅ Default value is clear (false/0/0.0)
5. ✅ Not dependent on complex logic

#### How to Add
```cpp
// In src/cell_struct.h, add to appropriate section:

// Boolean flags section:
bool my_new_flag{false};

// Integer properties section:
int my_new_count{0};

// Double properties section:
double my_new_value{0.0};
```

#### Example: Adding a Missing Field
```cpp
// If you see:
pCell->bIsProcessed()  // Field doesn't exist

// And it's clearly a simple boolean flag:
// 1. Add to cell_struct.h:
bool is_processed{false};

// 2. Use it:
cell.is_processed

// 3. Add comment in code:
// ADDED: is_processed field to Cell struct (line 142 in cell_struct.h)
```

#### When NOT to Add (Flag Instead)
```cpp
// DON'T add if:
// - Involves calculation/logic
// - Depends on other complex state
// - Purpose is unclear
// - Type is complex (vector, unique_ptr, etc.)

// Instead, FLAG:
// TODO AI-MIGRATION: Field 'pCell->SomeComplexMethod()' requires manual review
//   - Cannot determine appropriate field type/location
//   - May need refactoring rather than simple field addition
//   - See original at line XXX
```

---

## Loop Patterns

### Pattern A: Simple Grid-Wide Loop

```cpp
// BEFORE:
for (int nY = 0; nY < m_nYGridSize; nY++)
   for (int nX = 0; nX < m_nXGridSize; nX++)
      m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(0.0);

// AFTER:
auto& grid = m_pRasterGrid->Grid();
for (int nY = 0; nY < m_nYGridSize; nY++)
{
   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      grid(nX, nY).wave_height = 0.0;
   }
}
```

### Pattern B: Loop with Multiple Cell Accesses

```cpp
// BEFORE:
for (int nY = 0; nY < m_nYGridSize; nY++)
{
   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      double h = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();
      double a = m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle();
      double result = h * cos(a);
      m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth(result);
   }
}

// AFTER (cache cell reference):
auto& grid = m_pRasterGrid->Grid();
for (int nY = 0; nY < m_nYGridSize; nY++)
{
   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      Cell& cell = grid(nX, nY);
      double h = cell.wave_height;
      double a = cell.wave_angle;
      double result = h * cos(a);
      cell.sea_depth = result;
   }
}
```

### Pattern C: Range-Based Loop (When Coordinates Not Needed)

```cpp
// AFTER (modern - if you don't need x,y coordinates):
auto& grid = m_pRasterGrid->Grid();
for (Cell& cell : grid)
{
   cell.wave_height = 0.0;
}
```

---

## Flagging System

### When to Flag

Use `// TODO AI-MIGRATION:` comments for:

1. **Complex method calls** that involve logic
2. **Missing fields** that are complex or unclear
3. **Unusual access patterns** you don't understand
4. **Potential behavior changes** that need verification
5. **Compilation errors** you can't resolve
6. **Any uncertainty** about the correct migration

### How to Flag

```cpp
// TODO AI-MIGRATION: [Category] - Brief description
//   Context: Original code context
//   Issue: What's the problem
//   Location: File:Line reference if known
//   Recommendation: Suggested approach if any

// Example 1: Complex method
// TODO AI-MIGRATION: Complex Method - CalcSomething() involves non-trivial logic
//   Context: Was pCell->CalcComplexValue(param1, param2)
//   Issue: Method not present in Cell struct, logic unclear
//   Location: Original in cell.cpp:450
//   Recommendation: May need to become free function or grid method

// Example 2: Missing field
// TODO AI-MIGRATION: Missing Field - pCell->dGetSpecialValue()
//   Context: Used in erosion calculation
//   Issue: Field not in Cell struct, purpose unclear from context
//   Location: do_cliff_collapse.cpp:234
//   Recommendation: Need to determine if this is computed or stored

// Example 3: Uncertainty
// TODO AI-MIGRATION: Uncertain - Unusual pointer access pattern
//   Context: CGeomCell** ppCell = &m_pRasterGrid->m_Cell[nX][nY]
//   Issue: Double pointer dereference, unclear modern equivalent
//   Location: locate_coast.cpp:567
//   Recommendation: Manual review needed
```

### Flag Categories

Use these categories for consistency:
- `[Complex Method]` - Method with non-trivial logic
- `[Missing Field]` - Field not in Cell struct, complex to add
- `[Uncertain]` - Not sure about correct migration
- `[Potential Bug]` - Behavior may change, needs verification
- `[Performance]` - Access pattern may have performance implications
- `[Architecture]` - May require architectural changes

---

## Verification Steps

After each file migration:

### 1. Check Syntax
- Code compiles without errors
- No typos in field names
- Proper use of `.` vs `->` operators

### 2. Verify Field Names
- Cross-reference with `cell_struct.h`
- Check `member_mapping.txt` for old→new names
- Ensure snake_case, not camelCase

### 3. Check Reference Caching
- Grid reference cached when accessed multiple times
- Cell reference cached within loops
- No unnecessary repeated `Grid()` calls

### 4. Verify Access Patterns
- Old: `m_Cell[x][y]` → New: `Grid()(x, y)`
- Old: `->method()` → New: `.field` or `.method()`
- Nested: `pGetLandform()->` → `landform.`

### 5. Count TODOs
- List all `TODO AI-MIGRATION` comments added
- Provide summary of flagged issues by category

---

## Special Cases

### Case 1: pGetCell Returns Pointer

```cpp
// BEFORE:
CGeomCell* pCell = m_pRasterGrid->pGetCell(nX, nY);
pCell->dGetWaveHeight();

// AFTER (use reference):
Cell& cell = m_pRasterGrid->Grid()(nX, nY);
cell.wave_height;

// Only use pointer if truly needed:
Cell* pCell = &m_pRasterGrid->Grid()(nX, nY);
```

### Case 2: Const Access

```cpp
// BEFORE:
CGeomCell const* pCell = m_pRasterGrid->pGetCell(nX, nY);
double h = pCell->dGetWaveHeight();

// AFTER:
Cell const& cell = m_pRasterGrid->Grid()(nX, nY);
double h = cell.wave_height;
```

### Case 3: Conditional Access (Old Cell May Be Null)

```cpp
// BEFORE:
if (pCell)
   pCell->SetWaveHeight(5.0);

// AFTER (cells are never null in new system):
// Remove null check, grid always contains valid cells
cell.wave_height = 5.0;

// But FLAG if you're uncertain:
// TODO AI-MIGRATION: [Uncertain] - Null check removed
//   Context: Original had if (pCell) guard
//   Issue: Modern grid cells always valid, but verify this is correct
```

### Case 4: Grid Dimensions

```cpp
// BEFORE:
m_pRasterGrid->m_Cell[nX][nY].AppendLayers(n);

// AFTER:
m_pRasterGrid->Grid()(nX, nY).layers.resize(n);
```

### Case 5: Vector Operations on Layers

```cpp
// BEFORE:
int nLayers = pCell->nGetNumLayers();
pCell->AppendLayers(n);

// AFTER:
size_t nLayers = cell.layers.size();
cell.layers.resize(n);
```

---

## Example File Migration

### Before (simulation.cpp fragment)

```cpp
void CSimulation::UpdateWaves()
{
   for (int nY = 0; nY < m_nYGridSize; nY++)
   {
      for (int nX = 0; nX < m_nXGridSize; nX++)
      {
         CGeomCell* pCell = m_pRasterGrid->pGetCell(nX, nY);

         double dDeepHeight = pCell->dGetDeepWaterWaveHeight();
         double dDeepAngle = pCell->dGetDeepWaterWaveAngle();

         pCell->SetWaveHeight(dDeepHeight * 0.9);
         pCell->SetWaveAngle(dDeepAngle);

         if (pCell->bIsInContiguousSea())
         {
            pCell->pGetLandform()->SetAccumWaveEnergy(100.0);
         }
      }
   }
}
```

### After (Modern)

```cpp
void CSimulation::UpdateWaves()
{
   // Cache grid reference once
   auto& grid = m_pRasterGrid->Grid();

   for (int nY = 0; nY < m_nYGridSize; nY++)
   {
      for (int nX = 0; nX < m_nXGridSize; nX++)
      {
         // Cache cell reference for this iteration
         Cell& cell = grid(nX, nY);

         // Direct field access
         double dDeepHeight = cell.deep_water_wave_height;
         double dDeepAngle = cell.deep_water_wave_angle;

         cell.wave_height = dDeepHeight * 0.9;
         cell.wave_angle = dDeepAngle;

         if (cell.in_contiguous_sea)
         {
            cell.landform.accum_wave_energy = 100.0;
         }
      }
   }
}
```

---

## File Priority Order

Migrate in this order (from least to most complex):

1. **Initialization files** (low complexity)
   - `init_grid.cpp`
   - `read_input.cpp`

2. **Simple calculation files** (medium complexity)
   - `calc_curvature.cpp`
   - `calc_shadow_zones.cpp`

3. **Process files** (higher complexity)
   - `do_beach_sediment_movement.cpp`
   - `do_cliff_collapse.cpp`
   - `do_shore_platform_erosion.cpp`

4. **Complex files** (highest complexity, flag liberally)
   - `simulation.cpp` (large, many patterns)
   - `coast.cpp` (complex interactions)
   - `profile.cpp` (complex calculations)

---

## Quality Checklist

Before marking a file as complete:

- [ ] All `m_Cell[x][y]` replaced with `Grid()(x, y)`
- [ ] All getters/setters converted to direct field access
- [ ] Grid reference cached when accessed multiple times
- [ ] Cell reference cached in loops
- [ ] Field names verified against `cell_struct.h`
- [ ] No typos in field names (snake_case)
- [ ] Proper use of `.` vs `->` operators
- [ ] Complex cases flagged with `TODO AI-MIGRATION`
- [ ] Code compiles (or flags explain why not)
- [ ] No behavior changes introduced

---

## Communication Format

For each file migrated, report:

```markdown
## File: src/filename.cpp

**Status:** ✅ Complete | ⚠️ Partial | ❌ Blocked

**Changes:**
- Replaced X getter/setter calls
- Cached grid references in Y locations
- Cached cell references in Z loops
- Added N fields to cell_struct.h

**Flags Added:** M
- [Category] Line XXX: Brief description
- [Category] Line YYY: Brief description

**Compilation:** ✅ Passes | ❌ Errors (see flags)

**Confidence:** High | Medium | Low
```

---

## Example Workflow

For each file:

1. **Scan file** for old access patterns
2. **Read cell_struct.h** to know what fields exist
3. **Apply patterns** according to decision tree
4. **Cache references** for performance
5. **Flag complex cases** liberally
6. **Verify syntax** and field names
7. **Report results** with summary

---

## Common Mistakes to Avoid

❌ **Don't do this:**
```cpp
// Repeated Grid() calls:
m_pRasterGrid->Grid()(x, y).field1 = 1;
m_pRasterGrid->Grid()(x, y).field2 = 2;  // Wasteful!

// Wrong operator:
cell->wave_height  // cell is reference, use . not ->

// Typo in field name:
cell.waveHeight  // Should be snake_case: wave_height

// Add complex logic to Cell struct:
// (Flag instead!)
```

✅ **Do this:**
```cpp
// Cache reference:
Cell& cell = m_pRasterGrid->Grid()(x, y);
cell.field1 = 1;
cell.field2 = 2;

// Correct operator:
cell.wave_height

// Correct naming:
cell.wave_height

// Flag complex cases:
// TODO AI-MIGRATION: ...
```

---

## Success Criteria

Migration is successful when:

1. ✅ All simple field accesses converted
2. ✅ All nested component accesses converted
3. ✅ Code compiles with no errors
4. ✅ Performance optimizations applied (cached references)
5. ✅ Complex cases clearly flagged for manual review
6. ✅ No behavior changes introduced
7. ✅ Field names match cell_struct.h exactly

---

## Emergency: If Stuck

If you encounter something you can't handle:

1. **FLAG IT** - Add TODO AI-MIGRATION comment
2. **SKIP IT** - Move to next access pattern
3. **DOCUMENT** - Note in file report
4. **CONTINUE** - Don't block on one issue

Remember: **When in doubt, flag it out!**

Better to flag 100 cases for review than to introduce 1 bug.

---

## Final Note

You are modernizing critical infrastructure code. Prioritize:
1. **Correctness** over speed
2. **Clarity** over cleverness
3. **Flagging** over guessing

When uncertain, always choose to flag for manual review.

Good luck! 🚀
