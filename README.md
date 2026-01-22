# Blume-Capel-3D

A Monte Carlo simulation program for the Blume-Capel Model in three dimensions.

## Overview

This code simulates the **Blume-Capel Model in 3D**. It can also simulate the **Ising model** by setting `D = -1000` (or any large negative value), since for `D → -∞`, the Blume-Capel model converges to the Ising Model.

---

## Input Parameters

The program requires the following parameters to be entered via console:

### Lattice Size
- **Type:** `INT > 0`
- **Description:** Size of the 3D lattice

### Monte Carlo Steps
- **Type:** `INT > 0`
- **Recommendation for testing:** `10000`
- **Description:** Number of Monte Carlo steps to perform

### Number of Bins
- **Type:** `INT > 0` where `n_steps % N_Bin == 0`
- **Recommendation for testing:** `100`
- **Description:** Number of bins used in the Jackknife routine

### Model Selection

Choose between the **Blume-Capel** or **Ising** model:

#### Blume-Capel (BC)
- `β_c = 0.387721735`
- `D = 0.655`
- You can optionally change the value of `D`

#### Ising (I)
- `β_c = 0.22165463`

### High Temperature Mode

You can set a high temperature (e.g., `β = 0.1`). When enabled, the program will:
- Calculate the high-temperature expansion of the susceptibility for the BC model at `D = 0.641`
- Provide both the calculated value and a comparison to the measured value from the MC simulation

---

## Boundary Conditions

The program supports two types of boundary conditions for correlation calculations:

### 1. Periodic Boundary Conditions Only

**Selection:** Choose `'P'`

- All directions will have periodic boundary conditions
- Correlation function inside the bulk (bulk-bulk) will be calculated
- **Recommended for:**
  - High-temperature susceptibility expansion comparison
  - Ising model calculations
  - Binder Cumulant calculations

### 2. Open Boundary Conditions

**Selection:** Choose `'o'`

You can specify which directions have periodic boundary conditions. Correlations will only be calculated for surfaces in ±z direction.

**To disable periodic boundary conditions in ±z direction, enter:**
```
+x = 1
-x = 1
+y = 1
-y = 1
+z = 0
-z = 0
```

**Legend:**
- `1` = Periodic boundary conditions **ON**
- `0` = Periodic boundary conditions **OFF** (Open boundary conditions)

---

## External Field

**Prompt:** Do you want an external field `h ≠ 0`?

- Select **Yes** to enable
- Choose a value for the external field strength

---

## Additional Observables

**Prompt:** Calculate other observables? (Y/N)

Choose `'Y'` to calculate:
- Magnetization
- Susceptibility
- Binder Cumulant

> **Note:** This is enabled by default when calculating the high-temperature expansion.

---

## Important Notes

⚠️ **Correlation Functions:**
This program will always calculate correlations when the boundary conditions described above are chosen. These correlations may not generate meaningful results in certain scenarios, such as:
- When temperature is changed
- When magnetic field is enabled
- Other non-standard configurations

---

## Usage

1. Run the program
2. Follow the console prompts to input parameters
3. The simulation will execute and output results based on your configuration

---

## License

_Add your license information here_

## Contributing

_Add contributing guidelines here_

## Authors

_Add author information here_
