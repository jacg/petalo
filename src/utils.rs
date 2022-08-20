use std::error::Error;
use std::ops::{Bound, Range};

use units::{Ratio, todo::{Timef32, Lengthf32}};
use crate::{Point, BoundPair};
use crate::system_matrix::LOR;
use units::{mm, ns, ratio};

pub fn parse_range<T: std::str::FromStr>(s: &str) -> Result<Range<T>, <T as std::str::FromStr>::Err> {
    let v = s.split("..").collect::<Vec<_>>();
    if v.len() != 2 {
        panic!("Could not find '..' when parsing range.");
    }
    let x = v[0].parse()?;
    let y = v[1].parse()?;
    Ok(x..y)
}

pub fn parse_bounds<T: std::str::FromStr>(s: &str) -> Result<BoundPair<T>, <T as std::str::FromStr>::Err> {
    let v = s.split("..").collect::<Vec<_>>();
    if v.len() != 2 {
        panic!("Could not find '..' when parsing range.");
    }

    Ok((option_to_included(parse_if_not_empty(v[0])?),
        option_to_excluded(parse_if_not_empty(v[1])?)))

}

fn option_to_included<T>(x: Option<T>) -> Bound<T> {
    if let Some(y) = x { Bound::Included(y) }
    else               { Bound::Unbounded }
}

fn option_to_excluded<T>(x: Option<T>) -> Bound<T> {
    if let Some(y) = x { Bound::Excluded(y) }
    else               { Bound::Unbounded }
}

fn parse_if_not_empty<T: std::str::FromStr>(s: &str) -> Result<Option<T>, <T as std::str::FromStr>::Err> {
    Ok(if s.is_empty() { None }
       else            { Some(s.parse()?) })
}

#[allow(clippy::many_single_char_names)]
pub fn parse_triplet<T: std::str::FromStr>(s: &str) -> Result<(T,T,T), <T as std::str::FromStr>::Err> {
    let v = s.split(',').collect::<Vec<_>>();
    assert!(v.len() == 3);
    let x = v[0].parse()?;
    let y = v[1].parse()?;
    let z = v[2].parse()?;
    Ok((x, y, z))
}

pub fn parse_lor(s: &str) -> Result<LOR, Box<dyn Error>> {
    let n = s.split_whitespace().collect::<Vec<_>>();
    assert!(n.len() == 8);

    let t1 = ns(n[0].parse::<Timef32>()?);
    let t2 = ns(n[1].parse::<Timef32>()?);

    let x1 = mm(n[2].parse::<Lengthf32>()?);
    let y1 = mm(n[3].parse::<Lengthf32>()?);
    let z1 = mm(n[4].parse::<Lengthf32>()?);

    let x2 = mm(n[5].parse::<Lengthf32>()?);
    let y2 = mm(n[6].parse::<Lengthf32>()?);
    let z2 = mm(n[7].parse::<Lengthf32>()?);

    let p1 = Point::new(x1, y1, z1);
    let p2 = Point::new(x2, y2, z2);
    let lor = LOR::new(t1, t2, p1, p2, ratio(1.0));
    Ok(lor)
}

// Alias to disable structopt's type magic
pub type CutoffOption<T> = Option<T>;


pub fn parse_maybe_cutoff(s: &str) -> Result<CutoffOption<Ratio>, std::num::ParseFloatError> {
    Ok(if s == "no" { None } else { Some(units::ratio(s.parse()?)) })
}

/// Group numeric digits to facilitate reading long numbers
pub fn group_digits<F: std::fmt::Display>(n: F) -> String {
    use numsep::{separate, Locale};
    separate(n, Locale::English)
}


pub mod timing {

    use super::group_digits;
    use std::time::Instant;
    use std::io::Write;

    pub struct Progress {
        previous: Instant,
    }

    impl Progress {

        #[allow(clippy::new_without_default)]
        pub fn new() -> Self { Self { previous: Instant::now() } }

        /// Print message, append ellipsis, flush stdout, stay on same line, start timer.
        pub fn start(&mut self, message: &str) {
            print!("{message} ... ");
            std::io::stdout().flush().unwrap();
            self.start_timer();
        }

        /// Print message, go to next line, start timer
        pub fn startln(&mut self, message: &str) {
            self.start(message);
            println!();
            self.start_timer();
        }

        // Print time elapsed since last start or done
        pub fn done(&mut self) {
            println!("{} ms", group_digits(self.previous.elapsed().as_millis()));
            self.start_timer();
        }

        // Print message followed by time elapsed since last start or done
        pub fn done_with_message(&mut self, message: &str) {
            println!("{message}: {} ms",
                     group_digits(self.previous.elapsed().as_millis()));
            self.start_timer();
        }

        fn start_timer(&mut self) { self.previous = Instant::now() }
    }
}
