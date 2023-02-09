use units::{Length, PerLength, ratio_, mm, mm_, uom::ConstZero};

use geometry::Vector;

use crate::{
    LOR,
    config::mlem::Tof,
    fov::FOV,
    gauss::{make_gauss_option, Gaussian},
    system_matrix::{SystemMatrixRow, lor_fov_hit, FovHit, forward_project, back_project},
};

use super::{FoldState, Projector};

pub struct Siddon;

impl Projector for Siddon {

    fn project_one_lor<'i, 'g>(fold_state: FoldState<'i, 'g>, lor: &LOR) -> FoldState<'i, 'g> {
        project_one_lor(fold_state, lor)
    }

    fn buffers(fov: FOV) -> SystemMatrixRow {
        let [nx, ny, nz] = fov.n;
        let max_number_of_coupled_voxels_possible = nx + ny + nz - 2;
        SystemMatrixRow(Vec::with_capacity(max_number_of_coupled_voxels_possible))
    }

}

impl Siddon {
    // TODO  FOV and TOF should become construction-time arguments
    pub fn new_system_matrix_row(lor: &LOR, fov: &FOV, tof: Option<Tof>) -> SystemMatrixRow {
        let tof = make_gauss_option(tof);
        let mut system_matrix_row = Self::buffers(*fov);
        match lor_fov_hit(lor, *fov) {
            None => (),
            Some(FovHit {next_boundary, voxel_size, index, delta_index, remaining, tof_peak}) => {
                Self::update_smatrix_row(
                    &mut system_matrix_row,
                    next_boundary, voxel_size,
                    index, delta_index, remaining,
                    tof_peak, &tof
                );
            }
        }
        system_matrix_row
    }

    /// For a single LOR, place the weights and indices of the coupled voxels in
    /// `system_matrix_row` parameter. Using output parameters rather than return
    /// values, because this function is called in the inner loop, and allocating
    /// the vectors of results repeatedly, had a noticeable impact on performance.
    #[inline]
    #[allow(clippy::too_many_arguments)]
    pub fn update_smatrix_row(
        system_matrix_row: &mut SystemMatrixRow,
        mut next_boundary: Vector,
        voxel_size: Vector,
        mut index: i32,
        delta_index: [i32; 3],
        mut remaining: [i32; 3],
        tof_peak: Length,
        tof: &Option<Gaussian>
    ) {
        // Throw away previous LOR's values
        system_matrix_row.clear();

        // How far we have moved since entering the FOV
        let mut here = Length::ZERO;

        loop {
            // Which voxel boundary will be hit next, and its position
            let (dimension, boundary_position) = next_boundary.argmin();

            // The weight is the length of LOR in this voxel
            let mut weight = boundary_position - here;

            // If TOF enabled, adjust weight
            if let Some(gauss) = &tof {
                let g: PerLength = gauss.call(here - tof_peak);
                // TODO Normalization
                let completely_arbitrary_factor = 666.0;
                let g: f32 = ratio_(mm(completely_arbitrary_factor) * g);
                weight *= g;
            }

            // Store the index and weight of the voxel we have just crossed
            if weight > Length::ZERO {
                system_matrix_row.0.push(((index as usize), (mm_(weight))));
            }

            // Move along LOR until it leaves this voxel
            here = boundary_position;

            // Find the next boundary in this dimension
            next_boundary[dimension] += voxel_size[dimension];

            // Move index across the boundary we are crossing
            index += delta_index[dimension];
            remaining[dimension] -= 1;

            // If we have traversed the whole FOV, we're finished
            if remaining[dimension] == 0 { break; }
        }
    }

}

pub fn project_one_lor<'i, 'g>(state: FoldState<'i, 'g>, lor: &LOR) -> FoldState<'i, 'g> {
    let (mut backprojection, mut system_matrix_row, image, tof) = state;

    // Need to return the state from various match arms
    macro_rules! return_state { () => (return (backprojection, system_matrix_row, image, tof)); }

    // Analyse point where LOR hits FOV
    match lor_fov_hit(lor, image.fov) {

        // LOR missed FOV: nothing to be done
        None => return_state!(),

        // Data needed by `system_matrix_elements`
        Some(FovHit {next_boundary, voxel_size, index, delta_index, remaining, tof_peak}) => {

            // Find non-zero elements (voxels coupled to this LOR) and their values
            Siddon::update_smatrix_row(
                &mut system_matrix_row,
                next_boundary, voxel_size,
                index, delta_index, remaining,
                tof_peak, &tof
            );

            // Skip problematic LORs TODO: Is the cause more interesting than 'effiing floats'?
            for (i, _) in &system_matrix_row {
                if i >= backprojection.len() { return_state!(); }
            }

            // Forward projection of current image into this LOR
            let projection = forward_project(&system_matrix_row, image) * lor.additive_correction;

            // Backprojection of LOR onto image
            back_project(&mut backprojection, &system_matrix_row, ratio_(projection));
            return_state!();
        }
    }
}
