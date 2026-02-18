/// Parameters for Gasteiger electronegativity (a + bq + cq^2).
#[derive(Debug, Clone, Copy)]
pub struct GasteigerParams {
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

/// Supported hybridizations for different elements.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Hybridization {
    Sp3,
    Sp2,
    Sp,
    Default,
}

/// Helper function to get electronegativity parameters for an atom.
pub fn get_params(atomic_number: usize, hybridization: Hybridization) -> Option<GasteigerParams> {
    match (atomic_number, hybridization) {
        // Hydrogen
        (1, _) => Some(GasteigerParams { a: 7.17, b: 6.24, c: -0.56 }),
        
        // Carbon
        (6, Hybridization::Sp3) => Some(GasteigerParams { a: 7.98, b: 9.18, c: 1.88 }),
        (6, Hybridization::Sp2) => Some(GasteigerParams { a: 8.79, b: 9.32, c: 1.51 }),
        (6, Hybridization::Sp) => Some(GasteigerParams { a: 10.39, b: 9.45, c: 0.73 }),
        
        // Nitrogen
        (7, Hybridization::Sp3) => Some(GasteigerParams { a: 11.54, b: 10.82, c: 1.36 }),
        (7, Hybridization::Sp2) => Some(GasteigerParams { a: 12.87, b: 11.15, c: 0.85 }),
        (7, Hybridization::Sp) => Some(GasteigerParams { a: 15.68, b: 11.7, c: -0.27 }),
        
        // Oxygen
        (8, Hybridization::Sp3) => Some(GasteigerParams { a: 14.12, b: 12.92, c: 1.39 }),
        (8, Hybridization::Sp2) => Some(GasteigerParams { a: 17.07, b: 13.79, c: 0.47 }),
        
        // Fluorine
        (9, _) => Some(GasteigerParams { a: 14.66, b: 13.85, c: 2.31 }),
        
        // Chlorine
        (17, _) => Some(GasteigerParams { a: 10.18, b: 9.38, c: 1.13 }),
        
        // Bromine
        (35, _) => Some(GasteigerParams { a: 9.9, b: 8.29, c: 1.01 }),
        
        // Iodine
        (53, _) => Some(GasteigerParams { a: 8.85, b: 7.17, c: 0.99 }),

        // Phosphorus
        (15, _) => Some(GasteigerParams { a: 8.90, b: 8.12, c: 0.31 }),

        // Sulfur
        (16, Hybridization::Sp3) => Some(GasteigerParams { a: 10.14, b: 9.13, c: 1.38 }),
        (16, Hybridization::Sp2) => Some(GasteigerParams { a: 10.88, b: 9.47, c: 1.33 }),

        _ => None, // Fallback for unsupported elements/states
    }
}
