# gasteiger-rs

![Crates.io Version](https://img.shields.io/crates/v/gasteiger-rs)
![License](https://img.shields.io/badge/license-MIT%2FApache--2.0-blue)

A lightweight, zero-dependency Rust crate for calculating partial atomic charges using the **Gasteiger-Marsili PEOE (Partial Equalization of Orbital Electronegativity)** method.

This crate assigns partial charges based purely on molecular connectivity (topology), making it ideal for assigning initial charges before force field optimization (e.g., UFF) or for descriptors in chemoinformatics.

## Status: v0.9.0 (Beta)

This crate is currently in its **0.9.x** release phase. While the core algorithm (PEOE) and common elements (H, C, N, O, S, P, Halogens) are implemented and tested, results for complex or uncommon coordination environments should be verified. We are actively refining the hybridization detection logic and parameter tables.

## Features

- **Pure Rust & Zero Dependencies:** Built using only the standard library (`std`), ensuring fast compilation and easy integration.
- **Trait-Based Interface:** Works directly with your existing atom and bond structs via `GasteigerAtom` and `GasteigerBond` traits. No need to convert data structures.
- **Topology-Driven:** Calculates charges solely from the molecular graph. Does not require 3D coordinates.
- **Smart Hybridization Detection:** Automatically infers hybridization states ($sp^3, sp^2, sp$) based on connectivity and bond orders.
- **Robust Handling:**
    - **Ion Support:** Correctly distributes formal charges (e.g., $NH_4^+$, $CH_3COO^-$).
    - **Unknown Elements:** Safely handles unsupported elements by skipping charge transfer while preserving their formal charge.
    - **Sulfur/Phosphorus:** Special handling for hypervalent states.

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
gasteiger-rs = "0.9.0"
```

## Quick Start

Implement `GasteigerAtom` and `GasteigerBond` for your data structures and run the solver.

```rust
use gasteiger_rs::{GasteigerAtom, GasteigerBond, GasteigerSolver};

struct MyAtom {
    element: usize,
    formal_charge: f32,
}

impl GasteigerAtom for MyAtom {
    fn atomic_number(&self) -> usize { self.element }
    fn formal_charge(&self) -> f32 { self.formal_charge }
}

struct MyBond {
    pair: (usize, usize),
    order: f32,
}

impl GasteigerBond for MyBond {
    fn atom_indices(&self) -> (usize, usize) { self.pair }
    fn bond_order(&self) -> f32 { self.order }
}

fn main() {
    // Methane (CH4)
    let atoms = vec![
        MyAtom { element: 6, formal_charge: 0.0 }, // C
        MyAtom { element: 1, formal_charge: 0.0 }, // H
        MyAtom { element: 1, formal_charge: 0.0 }, // H
        MyAtom { element: 1, formal_charge: 0.0 }, // H
        MyAtom { element: 1, formal_charge: 0.0 }, // H
    ];
    let bonds = vec![
        MyBond { pair: (0, 1), order: 1.0 },
        MyBond { pair: (0, 2), order: 1.0 },
        MyBond { pair: (0, 3), order: 1.0 },
        MyBond { pair: (0, 4), order: 1.0 },
    ];

    let solver = GasteigerSolver::default();
    let charges = solver.compute_charges(&atoms, &bonds);

    println!("Partial Charges: {:?}", charges);
}
```

## Theory

The Partial Equalization of Orbital Electronegativity (PEOE) method iteratively transfers charge between bonded atoms to equalize their effective electronegativities.

1.  **Electronegativity ($\chi$):** Modeled as a quadratic function of charge $q$:
    $$ \chi = a + bq + cq^2 $$
2.  **Charge Transfer:** In each step, charge $dq$ moves from a less electronegative atom to a more electronegative one, damped by a factor $0.5^k$.
3.  **Convergence:** Typically converges within 6 iterations.

## License



MIT or Apache-2.0



## Author



**Forblaze Project**  

Website: [https://forblaze-works.com/](https://forblaze-works.com/)
