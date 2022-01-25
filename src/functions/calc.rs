type Point = (f64, f64, f64);

pub fn calc_angle(a: Point, b: Point, c: Point) -> f64 {
    // Form the two vectors
    let ba = [a.0 - b.0, a.1 - b.1, a.2 - b.2];
    let bc = [c.0 - b.0, c.1 - b.1, c.2 - b.2];

    // Calculate absolute values of vectors
    let abs_ba = ba.iter().fold(0.0, |acc, x| acc + (x * x)).sqrt();
    let abs_bc = bc.iter().fold(0.0, |acc, x| acc + (x * x)).sqrt();

    // Form dot product between vecs
    let dot = ba
        .iter()
        .zip(bc.iter())
        .fold(0.0, |acc, (a, b)| acc + (a * b));

    // Calculate angle from all ingredients
    (dot / (abs_ba * abs_bc)).acos().to_degrees()
}

pub fn calc_dihedral(a: Point, b: Point, c: Point, d: Point) -> f64 {
    // Form vectors
    let ba = [a.0 - b.0, a.1 - b.1, a.2 - b.2];
    let bc = [c.0 - b.0, c.1 - b.1, c.2 - b.2];
    let cb = [b.0 - c.0, b.1 - c.1, b.2 - c.2];
    let cd = [d.0 - c.0, d.1 - c.1, d.2 - c.2];

    // Form two normal vectors via cross products
    let n1 = [
        ba[1] * bc[2] - ba[2] * bc[1],
        ba[2] * bc[0] - ba[0] * bc[2],
        ba[0] * bc[1] - ba[1] * bc[0],
    ];
    let n2 = [
        cb[1] * cd[2] - cb[2] * cd[1],
        cb[2] * cd[0] - cb[0] * cd[2],
        cb[0] * cd[1] - cb[1] * cd[0],
    ];

    // calculate abs of vecs
    let abs_n1 = n1.iter().fold(0.0, |acc, x| acc + (x * x)).sqrt();
    let abs_n2 = n2.iter().fold(0.0, |acc, x| acc + (x * x)).sqrt();

    let dot = n1
        .iter()
        .zip(n2.iter())
        .fold(0.0, |acc, (a, b)| acc + (a * b));
    (dot / (abs_n1 * abs_n2)).acos().to_degrees()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn angle_test() {
        let a = (1.0, 0.0, 0.0);
        let b = (0.0, 1.0, 0.0);
        let c = (0.0, 0.0, 1.0);

        assert!((60.0 - calc_angle(a, b, c)).abs() < 0.0001)
    }

    #[test]
    fn dihedral_test() {
        let a = (1.0, 0.0, 0.0);
        let b = (0.0, 0.0, 0.0);
        let c = (0.0, 0.0, 1.0);
        let d = (1.0, 1.0, 1.0);

        println!("{}", calc_dihedral(a, b, c, d));
        assert!((45.0 - calc_dihedral(a, b, c, d)).abs() < 0.0001);
    }
}
