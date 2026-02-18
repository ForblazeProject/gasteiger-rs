pub mod traits;
pub mod parameters;
pub mod solver;

pub use traits::{GasteigerAtom, GasteigerBond};
pub use solver::GasteigerSolver;

#[cfg(test)]
mod tests {
    use super::*;

    struct MockAtom {
        name: &'static str,
        element: usize,
        formal_charge: f32,
    }

    impl GasteigerAtom for MockAtom {
        fn atomic_number(&self) -> usize { self.element }
        fn formal_charge(&self) -> f32 { self.formal_charge }
    }

    struct MockBond {
        pair: (usize, usize),
        order: f32,
    }

    impl GasteigerBond for MockBond {
        fn atom_indices(&self) -> (usize, usize) { self.pair }
        fn bond_order(&self) -> f32 { self.order }
    }

    fn print_charges(atoms: &[MockAtom], charges: &[f64]) {
        println!("{:<5} {:<10} {:<15}", "Idx", "Element", "Partial Charge");
        for (i, (atom, &charge)) in atoms.iter().zip(charges.iter()).enumerate() {
            println!("{:<5} {:<10} {:<15.6}", i, atom.name, charge);
        }
        let total: f64 = charges.iter().sum();
        println!("Total Charge: {:.6}", total);
    }

    #[test]
    fn test_triple_bonds() {
        let solver = GasteigerSolver::default();

        // 1. Acetylene (HC#CH)
        let atoms_ace = vec![
            MockAtom { name: "C1", element: 6, formal_charge: 0.0 }, // Sp
            MockAtom { name: "C2", element: 6, formal_charge: 0.0 }, // Sp
            MockAtom { name: "H1", element: 1, formal_charge: 0.0 },
            MockAtom { name: "H2", element: 1, formal_charge: 0.0 },
        ];
        let bonds_ace = vec![
            MockBond { pair: (0, 1), order: 3.0 }, // Triple bond
            MockBond { pair: (0, 2), order: 1.0 },
            MockBond { pair: (1, 3), order: 1.0 },
        ];
        let charges_ace = solver.compute_charges(&atoms_ace, &bonds_ace);
        println!("\n--- Acetylene (HC#CH) ---");
        print_charges(&atoms_ace, &charges_ace);
        assert!(charges_ace.iter().sum::<f64>().abs() < 1e-6);

        // 2. Acetonitrile (CH3-C#N)
        let atoms_nit = vec![
            MockAtom { name: "C_me", element: 6, formal_charge: 0.0 }, // Sp3
            MockAtom { name: "C_cn", element: 6, formal_charge: 0.0 }, // Sp
            MockAtom { name: "N",    element: 7, formal_charge: 0.0 }, // Sp
            MockAtom { name: "H1",   element: 1, formal_charge: 0.0 },
            MockAtom { name: "H2",   element: 1, formal_charge: 0.0 },
            MockAtom { name: "H3",   element: 1, formal_charge: 0.0 },
        ];
        let bonds_nit = vec![
            MockBond { pair: (0, 1), order: 1.0 },
            MockBond { pair: (1, 2), order: 3.0 }, // C#N
            MockBond { pair: (0, 3), order: 1.0 },
            MockBond { pair: (0, 4), order: 1.0 },
            MockBond { pair: (0, 5), order: 1.0 },
        ];
        let charges_nit = solver.compute_charges(&atoms_nit, &bonds_nit);
        println!("\n--- Acetonitrile (CH3-CN) ---");
        print_charges(&atoms_nit, &charges_nit);
        assert!(charges_nit.iter().sum::<f64>().abs() < 1e-6);
    }

    #[test]
    fn test_ethene_and_benzene() {
        let solver = GasteigerSolver::default();

        let atoms_ethene = vec![
            MockAtom { name: "C1", element: 6, formal_charge: 0.0 },
            MockAtom { name: "C2", element: 6, formal_charge: 0.0 },
            MockAtom { name: "H1", element: 1, formal_charge: 0.0 },
            MockAtom { name: "H2", element: 1, formal_charge: 0.0 },
            MockAtom { name: "H3", element: 1, formal_charge: 0.0 },
            MockAtom { name: "H4", element: 1, formal_charge: 0.0 },
        ];
        let bonds_ethene = vec![
            MockBond { pair: (0, 1), order: 2.0 },
            MockBond { pair: (0, 2), order: 1.0 },
            MockBond { pair: (0, 3), order: 1.0 },
            MockBond { pair: (1, 4), order: 1.0 },
            MockBond { pair: (1, 5), order: 1.0 },
        ];
        let charges_ethene = solver.compute_charges(&atoms_ethene, &bonds_ethene);
        println!("\n--- Ethene (CH2=CH2) ---");
        print_charges(&atoms_ethene, &charges_ethene);
        assert!(charges_ethene.iter().sum::<f64>().abs() < 1e-6);

        let mut atoms_benzene = Vec::new();
        for _ in 0..6 { atoms_benzene.push(MockAtom { name: "C", element: 6, formal_charge: 0.0 }); }
        for _ in 0..6 { atoms_benzene.push(MockAtom { name: "H", element: 1, formal_charge: 0.0 }); }
        
        let mut bonds_benzene = Vec::new();
        for i in 0..6 {
            bonds_benzene.push(MockBond { pair: (i, (i + 1) % 6), order: 1.5 });
            bonds_benzene.push(MockBond { pair: (i, i + 6), order: 1.0 });
        }
        let charges_benzene = solver.compute_charges(&atoms_benzene, &bonds_benzene);
        println!("\n--- Benzene (C6H6) ---");
        print_charges(&atoms_benzene, &charges_benzene);
        assert!(charges_benzene.iter().sum::<f64>().abs() < 1e-6);
    }

    #[test]
    fn test_sulfur_compounds() {
        let solver = GasteigerSolver::default();
        let atoms = vec![
            MockAtom { name: "S", element: 16, formal_charge: 0.0 },
            MockAtom { name: "C1", element: 6, formal_charge: 0.0 },
            MockAtom { name: "C2", element: 6, formal_charge: 0.0 },
        ];
        let bonds = vec![
            MockBond { pair: (0, 1), order: 1.0 },
            MockBond { pair: (0, 2), order: 1.0 },
        ];
        let charges = solver.compute_charges(&atoms, &bonds);
        println!("\n--- Dimethyl Sulfide (CH3-S-CH3) ---");
        print_charges(&atoms, &charges);

        let atoms_thio = vec![
            MockAtom { name: "S", element: 16, formal_charge: 0.0 },
            MockAtom { name: "C", element: 6, formal_charge: 0.0 },
        ];
        let bonds_thio = vec![
            MockBond { pair: (0, 1), order: 2.0 },
        ];
        let charges_thio = solver.compute_charges(&atoms_thio, &bonds_thio);
        println!("\n--- Thioformaldehyde (CH2=S) ---");
        print_charges(&atoms_thio, &charges_thio);
        assert!(charges_thio[0] < 0.0);
    }

    #[test]
    fn test_unknown_element_handling() {
        let atoms = vec![
            MockAtom { name: "C", element: 6, formal_charge: 0.0 },
            MockAtom { name: "H1", element: 1, formal_charge: 0.0 },
            MockAtom { name: "Pd", element: 46, formal_charge: 0.0 },
        ];
        let bonds = vec![
            MockBond { pair: (0, 1), order: 1.0 },
            MockBond { pair: (0, 2), order: 1.0 },
        ];
        let solver = GasteigerSolver::default();
        let charges = solver.compute_charges(&atoms, &bonds);
        println!("\n--- Unknown Element Test (C-H + C-Pd) ---");
        print_charges(&atoms, &charges);
        assert!(charges[2].abs() < 1e-10);
    }

    #[test]
    fn test_methane_charges() {
        let atoms = vec![
            MockAtom { name: "C", element: 6, formal_charge: 0.0 },
            MockAtom { name: "H1", element: 1, formal_charge: 0.0 },
            MockAtom { name: "H2", element: 1, formal_charge: 0.0 },
            MockAtom { name: "H3", element: 1, formal_charge: 0.0 },
            MockAtom { name: "H4", element: 1, formal_charge: 0.0 },
        ];
        let bonds = vec![
            MockBond { pair: (0, 1), order: 1.0 },
            MockBond { pair: (0, 2), order: 1.0 },
            MockBond { pair: (0, 3), order: 1.0 },
            MockBond { pair: (0, 4), order: 1.0 },
        ];
        let solver = GasteigerSolver::default();
        let charges = solver.compute_charges(&atoms, &bonds);
        println!("\n--- Methane (CH4) ---");
        print_charges(&atoms, &charges);
        assert!(charges[0] < 0.0);
    }

    #[test]
    fn test_water_charges() {
        let atoms = vec![
            MockAtom { name: "O", element: 8, formal_charge: 0.0 },
            MockAtom { name: "H1", element: 1, formal_charge: 0.0 },
            MockAtom { name: "H2", element: 1, formal_charge: 0.0 },
        ];
        let bonds = vec![
            MockBond { pair: (0, 1), order: 1.0 },
            MockBond { pair: (0, 2), order: 1.0 },
        ];
        let solver = GasteigerSolver::default();
        let charges = solver.compute_charges(&atoms, &bonds);
        println!("\n--- Water (H2O) ---");
        print_charges(&atoms, &charges);
        assert!(charges[0] < 0.0);
    }

    #[test]
    fn test_multi_molecule_system() {
        let atoms = vec![
            MockAtom { name: "C", element: 6, formal_charge: 0.0 },
            MockAtom { name: "H1", element: 1, formal_charge: 0.0 },
            MockAtom { name: "H2", element: 1, formal_charge: 0.0 },
            MockAtom { name: "H3", element: 1, formal_charge: 0.0 },
            MockAtom { name: "H4", element: 1, formal_charge: 0.0 },
            MockAtom { name: "O", element: 8, formal_charge: 0.0 },
            MockAtom { name: "Hw1", element: 1, formal_charge: 0.0 },
            MockAtom { name: "Hw2", element: 1, formal_charge: 0.0 },
        ];
        let bonds = vec![
            MockBond { pair: (0, 1), order: 1.0 },
            MockBond { pair: (0, 2), order: 1.0 },
            MockBond { pair: (0, 3), order: 1.0 },
            MockBond { pair: (0, 4), order: 1.0 },
            MockBond { pair: (5, 6), order: 1.0 },
            MockBond { pair: (5, 7), order: 1.0 },
        ];
        let solver = GasteigerSolver::default();
        let charges = solver.compute_charges(&atoms, &bonds);
        println!("\n--- Multi-molecule (CH4 + H2O) ---");
        print_charges(&atoms, &charges);
        let methane_sum: f64 = charges[0..5].iter().sum();
        let water_sum: f64 = charges[5..8].iter().sum();
        assert!(methane_sum.abs() < 1e-6);
        assert!(water_sum.abs() < 1e-6);
    }

    #[test]
    fn test_ion_charge_conservation() {
        let atoms = vec![
            MockAtom { name: "N+", element: 7, formal_charge: 1.0 },
            MockAtom { name: "H1", element: 1, formal_charge: 0.0 },
            MockAtom { name: "H2", element: 1, formal_charge: 0.0 },
            MockAtom { name: "H3", element: 1, formal_charge: 0.0 },
            MockAtom { name: "H4", element: 1, formal_charge: 0.0 },
        ];
        let bonds = vec![
            MockBond { pair: (0, 1), order: 1.0 },
            MockBond { pair: (0, 2), order: 1.0 },
            MockBond { pair: (0, 3), order: 1.0 },
            MockBond { pair: (0, 4), order: 1.0 },
        ];
        let solver = GasteigerSolver::default();
        let charges = solver.compute_charges(&atoms, &bonds);
        println!("\n--- Ammonium Ion (NH4+) ---");
        print_charges(&atoms, &charges);
        let total_charge: f64 = charges.iter().sum();
        assert!((total_charge - 1.0).abs() < 1e-6);
    }
}
