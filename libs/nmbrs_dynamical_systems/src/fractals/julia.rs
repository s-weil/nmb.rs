use nalgebra::{Complex, Normed};

fn julia(c: Complex<f64>, z: Complex<f64>) -> Complex<f64> {
    z * z + c
}

const MAX_ITERATION: u64 = 255;

// fn julia_orbit(c: Complex<f64>, z: Complex<f64>, n: u64) -> Option<u64> {
//     if z.norm() > 2.0 || n >= MAX_ITERATION {
//         return Some(n);
//     }

//     julia_orbit(c, z, n + 1)
// }

fn julia_orbit(c: Complex<f64>, z: Complex<f64>) -> Option<u64> {
    let mut z = z;
    for i in 0..255 {
        if z.norm() > 2.0 {
            return Some(i);
        }
        z = julia(c, z);
    }

    None
}

#[cfg(test)]
mod tests {
    use image::{ImageBuffer, Rgb};
    use nalgebra::Complex;

    use super::{julia_orbit, MAX_ITERATION};

    #[test]
    fn plot() {
        let width: u32 = 800;
        let height = 600;

        let scale_x = 3. / width as f64;
        let scale_y = 3. / height as f64;

        let c = Complex::new(-0.8, 0.156);

        let mut img_buffer: ImageBuffer<_, Vec<_>> = ImageBuffer::new(width, height);

        for (x, y, pixel) in img_buffer.enumerate_pixels_mut() {
            let z = Complex::new(x as f64 * scale_x - 1.5, y as f64 * scale_y - 1.5);
            let n_iter = julia_orbit(c, z).unwrap_or(MAX_ITERATION) as u8;

            *pixel = Rgb([n_iter, n_iter, n_iter]);
        }

        img_buffer.save("julia.png").unwrap();
    }
}
