use std::error::Error;
use std::ops::Range;

use crate::types::{Time, Length};
use crate::weights::{LOR, Point};

pub fn parse_range<T: std::str::FromStr>(s: &str) -> Result<Range<T>, <T as std::str::FromStr>::Err> {
    let v = s.split("..").collect::<Vec<_>>();
    assert!(v.len() == 2);
    let x = v[0].parse()?;
    let y = v[1].parse()?;
    Ok(x..y)
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

    let t1 = n[0].parse::<Time>()?;
    let t2 = n[1].parse::<Time>()?;

    let x1 = n[2].parse::<Length>()?;
    let y1 = n[3].parse::<Length>()?;
    let z1 = n[4].parse::<Length>()?;

    let x2 = n[5].parse::<Length>()?;
    let y2 = n[6].parse::<Length>()?;
    let z2 = n[7].parse::<Length>()?;

    let p1 = Point::new(x1, y1, z1);
    let p2 = Point::new(x2, y2, z2);
    let lor = LOR::new(t1, t2, p1, p2);
    Ok(lor)
}
