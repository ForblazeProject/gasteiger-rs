/// Trait for atoms that need Gasteiger partial charges.
pub trait GasteigerAtom {
    /// Atomic number (e.g., H=1, C=6).
    fn atomic_number(&self) -> usize;
    /// Formal charge of the atom (default is 0.0).
    fn formal_charge(&self) -> f32 {
        0.0
    }
}

/// Trait for bonds between atoms.
pub trait GasteigerBond {
    /// Indices of the two atoms connected by this bond.
    fn atom_indices(&self) -> (usize, usize);
    /// Bond order (1.0 for single, 2.0 for double, 3.0 for triple, 1.5 for aromatic).
    fn bond_order(&self) -> f32;
}
