use crate::traits::{GasteigerAtom, GasteigerBond};
use crate::parameters::{Hybridization, get_params, GasteigerParams};

pub struct GasteigerSolver {
    pub iterations: usize,
    pub damping: f64,
}

impl Default for GasteigerSolver {
    fn default() -> Self {
        Self {
            iterations: 6,
            damping: 0.5,
        }
    }
}

impl GasteigerSolver {
    pub fn compute_charges<A, B>(&self, atoms: &[A], bonds: &[B]) -> Vec<f64>
    where
        A: GasteigerAtom,
        B: GasteigerBond,
    {
        let n_atoms = atoms.len();
        let mut charges = vec![0.0; n_atoms];
        
        for (i, atom) in atoms.iter().enumerate() {
            charges[i] = atom.formal_charge() as f64;
        }

        let mut atom_params: Vec<Option<GasteigerParams>> = Vec::with_capacity(n_atoms);
        for i in 0..n_atoms {
            let hybrid = self.guess_hybridization(i, atoms, bonds);
            let params = get_params(atoms[i].atomic_number(), hybrid)
                .or_else(|| get_params(atoms[i].atomic_number(), Hybridization::Sp3))
                .or_else(|| get_params(atoms[i].atomic_number(), Hybridization::Default));
            atom_params.push(params);
        }

        let mut current_damping = 1.0;
        for _ in 0..self.iterations {
            let mut delta_charges = vec![0.0; n_atoms];

            for bond in bonds {
                let (i, j) = bond.atom_indices();
                if i >= n_atoms || j >= n_atoms { continue; }

                if let (Some(p_i), Some(p_j)) = (&atom_params[i], &atom_params[j]) {
                    let chi_i = self.calculate_electronegativity(p_i, charges[i]);
                    let chi_j = self.calculate_electronegativity(p_j, charges[j]);

                    let chi_plus_i = self.calculate_electronegativity(p_i, 1.0);
                    let chi_plus_j = self.calculate_electronegativity(p_j, 1.0);

                    if chi_j > chi_i {
                        let dq = current_damping * (chi_j - chi_i) / chi_plus_i;
                        delta_charges[i] += dq;
                        delta_charges[j] -= dq;
                    } else if chi_i > chi_j {
                        let dq = current_damping * (chi_i - chi_j) / chi_plus_j;
                        delta_charges[j] += dq;
                        delta_charges[i] -= dq;
                    }
                }
            }

            for i in 0..n_atoms {
                charges[i] += delta_charges[i];
            }
            current_damping *= self.damping;
        }

        charges
    }

    fn calculate_electronegativity(&self, p: &GasteigerParams, q: f64) -> f64 {
        p.a + p.b * q + p.c * q * q
    }

    fn guess_hybridization<A, B>(&self, atom_idx: usize, atoms: &[A], bonds: &[B]) -> Hybridization
    where
        A: GasteigerAtom,
        B: GasteigerBond,
    {
        let atomic_number = atoms[atom_idx].atomic_number();
        let mut neighbor_count = 0;

        for bond in bonds {
            let (i, j) = bond.atom_indices();
            if i == atom_idx || j == atom_idx {
                neighbor_count += 1;
            }
        }

        match atomic_number {
            6 => { // Carbon
                if neighbor_count >= 4 { Hybridization::Sp3 }
                else if neighbor_count == 3 { Hybridization::Sp2 }
                else if neighbor_count <= 2 { Hybridization::Sp }
                else { Hybridization::Sp3 }
            }
            7 => { // Nitrogen
                if neighbor_count >= 3 { Hybridization::Sp3 }
                else if neighbor_count == 2 { Hybridization::Sp2 }
                else { Hybridization::Sp }
            }
            8 => { // Oxygen
                if neighbor_count >= 2 { Hybridization::Sp3 }
                else { Hybridization::Sp2 }
            }
            15 => { // Phosphorus
                Hybridization::Sp3
            }
            16 => { // Sulfur
                if neighbor_count >= 2 { Hybridization::Sp3 }
                else { Hybridization::Sp2 }
            }
            _ => Hybridization::Default,
        }
    }
}