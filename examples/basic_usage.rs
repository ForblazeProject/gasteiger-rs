use gasteiger_rs::{GasteigerAtom, GasteigerBond, GasteigerSolver};

#[derive(Debug)]
struct MyAtom {
    name: &'static str,
    element: usize,
    formal_charge: f32,
    partial_charge: f64,
}

impl GasteigerAtom for MyAtom {
    fn atomic_number(&self) -> usize {
        self.element
    }
    fn formal_charge(&self) -> f32 {
        self.formal_charge
    }
}

struct MyBond {
    pair: (usize, usize),
    order: f32,
}

impl GasteigerBond for MyBond {
    fn atom_indices(&self) -> (usize, usize) {
        self.pair
    }
    fn bond_order(&self) -> f32 {
        self.order
    }
}

fn main() {
    // 1. Define Atoms (Methane: CH4)
    let mut atoms = vec![
        MyAtom { name: "C",  element: 6, formal_charge: 0.0, partial_charge: 0.0 },
        MyAtom { name: "H1", element: 1, formal_charge: 0.0, partial_charge: 0.0 },
        MyAtom { name: "H2", element: 1, formal_charge: 0.0, partial_charge: 0.0 },
        MyAtom { name: "H3", element: 1, formal_charge: 0.0, partial_charge: 0.0 },
        MyAtom { name: "H4", element: 1, formal_charge: 0.0, partial_charge: 0.0 },
    ];

    // 2. Define Bonds
    let bonds = vec![
        MyBond { pair: (0, 1), order: 1.0 },
        MyBond { pair: (0, 2), order: 1.0 },
        MyBond { pair: (0, 3), order: 1.0 },
        MyBond { pair: (0, 4), order: 1.0 },
    ];

    // 3. Compute Charges
    let solver = GasteigerSolver::default();
    let charges = solver.compute_charges(&atoms, &bonds);

    // 4. Update Atoms with Results
    for (i, charge) in charges.iter().enumerate() {
        atoms[i].partial_charge = *charge;
        println!("Atom {:<2} ({:<2}): {:.6}", i, atoms[i].name, charge);
    }

    // Verify Total Charge
    let total_charge: f64 = charges.iter().sum();
    println!("Total Charge: {:.6}", total_charge);
}
