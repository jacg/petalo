use uom::fmt::DisplayStyle::Abbreviation;
use uom::si::f32::{Length, Time, Velocity};
use uom::si::{length  ::{kilometer  as km,
                         meter      as  m,
                         millimeter as mm},
              velocity::{meter_per_second},
              time    ::{second     as  s}};


fn main() {
    let c: Velocity = Velocity::new::<meter_per_second>(299_792_458.0);

    let l1: Length = Length::new::<km> (3.2);
    let l2         = Length::new::< m>(40.2);
    let l3         = Length::new::<mm>(99.2);
    let t          = Time::new::<s>(13.3);
    let v1: Velocity = (l1 + l2) / t;
    let v2           =       l2  / t;
    //let type_error = l1 + t;

    // Reusable format arguments
    let fmm = Length::format_args(mm, Abbreviation);
    let  fm = Length::format_args( m, Abbreviation);
    let fkm = Length::format_args(km, Abbreviation);

    println!("{} = {}",  fm.with(l1), fkm.with(l1));
    println!("{} = {}", fmm.with(l3),  fm.with(l3));
    println!("{}", v1.into_format_args(meter_per_second, Abbreviation));
    println!("{:?}", v2 / c);

}
