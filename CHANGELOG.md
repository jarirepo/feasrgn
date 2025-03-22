# Changelog

## [1.2.0][master] - 2016-12-08
### Changed
- Improved algorithm and visualization using the OOP capabilities.
- Renamed class to `inequality`
- Removed the auto-indexing
- Replaced `lt()` and `gt()` with `isNeg()` and `isPos()` respectively
- The signature of the inequality is obtained by `getSignature()`
- Default `lhsName` and `rhsName` are set by `setDefaultNames()`
- The logics related to the inequality sign is now handled by the enumeration class `ineqsign`
- Added enumeration class `ineqsign` for the inequality sign
- The inequalities now need to be defined in a cell array which allows feasrgn to automatically number the inequalities according to their order of occurence in the cell array
- Improved the scanning algorithm further by utilizing node connectivity information to ensure that all intersections will be included in the final generated boundary
- Colinear boundary segments are automatically removed
- The cell-array with legend string (legStr) is created in feasrgnplot to keep this class clean from the visualization
- Improved validation of the input
- New method isPtInside to test if a query point q is inside the feasible region boundary
- New method generateInnerPts to generate a 2-dimensional set of points inside the feasible region boundary B
- Class `feasrgnplot` for visualization of the constraints functions and the feasible region. Most style properties can be dynamically specified.

## [1.1.0] - 2016-11-02
### Changed
- Rewrote portions of the scanning algorithm which resolved some issues with missing intersections along the generated feasible region boundary and some other issues. The updated algorithm is more robust and also faster.

## [1.0.0] - 2016-11-01
### Added
- Initial algorithm

[master]: https://se.mathworks.com/matlabcentral/fileexchange/60049-jarirepo-feasrgn?s_tid=FX_rc1_behav
